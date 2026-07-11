/* MYCN isoform-stack viewer — custom SVG (d3 v7), superset of the igv viewer.
   Transcript-first isoform stack (fixed order, accordion ORFs, specificity ranking) PLUS
   chromosome ideogram, kb ruler, rainbow reference strip, zoom/pan, base-level nucleotides +
   codons, dual browse lens, search, collapsible drawers, now-showing marker, SVG + PNG export.
   All structure/ranking derived from mycn_orfs.meta.json; base letters from mycn_locus_seq.json
   (index = pos - 15,938,000). Nothing hand-typed; no coordinate re-derivation. Palette + card
   format match the igv viewer. */

(function () {
  "use strict";

  // ---- palette (loaded from META.palette at boot; these are v3.10 fallbacks) ----
  var PAL = { nmyc:"#4A6FA5", dmycn:"#D98A3D", mycnot:"#C0504D", musep:"#B0568C",
              ri_line:"#2F7A57", ri_bg:"#C3E0D4", pred_ext:"#8FA9C9", utr:"#C7CFCC",
              term:"#8A9791", gold:"#C9A227" };
  var BACKBONE = "#b9c2c0";
  var FONT = "Inter, system-ui, -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif";
  var MONO = "ui-monospace, 'SF Mono', 'Roboto Mono', Menlo, Consolas, monospace";
  // IGV-style nucleotide colours (A green, C blue, G orange, T red)
  var BASECOL = { A:"#33A02C", C:"#1F78B4", G:"#FF7F00", T:"#E31A1C", U:"#E31A1C", N:"#B0B0B0" };

  var TX_ORDER = ["201","202","203","204","205","206","207","XM"];
  var NAME = { "1":"N-Myc", "2":"ΔMYCN", "9":"MYCNOT", "10":"MUSEP" };
  // by-ORF browse groups (membership derived from cls; never hand-typed)
  var GROUPS = [
    { key:"canonical-CDS",     title:"Core proteins" },
    { key:"known-uORF",        title:"Known functional uORFs" },
    { key:"N-Myc-frame-ext",   title:"Predicted N-terminal extensions (206)" },
    { key:"N-Myc-frame-dsATG", title:"In-frame downstream-ATG ORFs" }
  ];
  var NOVEL = [ { key:"novel-uORF", title:"novel-uORF" }, { key:"novel-other", title:"novel-other" } ];
  var LEGEND = [];   // derived in buildLegendData() from PAL

  // ---- geometry ----
  var GUTTER = 200, PADR = 18;   // wider label lane so ORF-row columns (swatch/id+name/aa·carrier) never collide
  var TXH = 50, ORFH = 23, GRPH = 20, SHH = 27, EXPAD = 8;
  var THICK = 15, THIN = 6, RIH = 9, ORF_THICK = 11;
  var DMIN = 15938000, DMAX = 15947200;       // locus extent (matches staged sequence)
  var CHR2LEN = 242193529;                     // chr2 length (fallback; real value from CYTO.length)
  var LETTER_PXPB = 7, CELL_PXPB = 0.5;       // zoom thresholds (px per base)
  var HDR_H = 56, IDEO_H = 56;   // room below the bar for the band label (p24.3) — was 46 (clipped its descender)
  // Giemsa greyscale for cytoband stains (UCSC convention)
  var STAIN = { gneg:"#ffffff", gpos25:"#c9cdd1", gpos50:"#9aa1a7", gpos75:"#6f767c",
                gpos100:"#4a5157", gvar:"#c2c9cf", stalk:"#9aa7b0", acen:"#b06a80" };

  var META = null, LOCSEQ = null, CYTO = null;
  var x0 = null, x = null, zt = null, zoom = null, rulerBrush = null, dragPrevY = null;
  var svg = null, hdr = null, ideo = null, gMain = null;
  var W = 900, H = 400;
  var ROWS = [];
  var lens = "tx", selectedOrf = null, selectedTx = null, novelOpen = false, railHasContent = false, rightOpen = false;
  var viewG = [DMIN, DMAX];
  var state = {};
  var CANON_INTRONS = [], RETAINED = {}, CDS_PHASE = {}, GPOS = {};

  function $(id){ return document.getElementById(id); }
  function d3s(id){ return d3.select("#"+id); }
  function esc(s){ return String(s).replace(/[&<>"]/g, function(c){
    return {"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;"}[c]; }); }
  function fmt(n){ return Math.round(n).toLocaleString("en-US"); }
  function clamp(v,a,b){ return v<a?a:(v>b?b:v); }
  function toast(msg, isErr){ var el=$("loadmsg"); if(!el) return;
    if(!msg){ el.classList.add("hide"); return; }
    el.textContent=msg; el.classList.remove("hide"); el.classList.toggle("err", !!isErr); }

  // ================= boot =================
  window.addEventListener("DOMContentLoaded", function () {
    Promise.all([
      fetch("mycn_orfs.meta.json").then(function(r){ if(!r.ok) throw new Error("meta.json HTTP "+r.status); return r.json(); }),
      fetch("mycn_locus_seq.json").then(function(r){ if(!r.ok) throw new Error("locus_seq HTTP "+r.status); return r.json(); }),
      fetch("chr2_cytoband.json").then(function(r){ if(!r.ok) throw new Error("cytoband HTTP "+r.status); return r.json(); })
    ]).then(function(res){
      META = res[0]; LOCSEQ = res[1]; CYTO = res[2];
      if (!verifyBaseMap()) throw new Error("base-map ATG spot-check FAILED at chr2:15,942,065");
      init(); toast(null);
    }).catch(function(e){ console.error(e); toast("Failed to load MYCN data: "+e.message, true); });
  });

  // base-map integrity: canonical ATG at genomic 15,942,065 must read ATG (index = pos - DMIN)
  function verifyBaseMap(){
    if (!LOCSEQ || LOCSEQ.start !== DMIN) return false;
    var i = 15942065 - LOCSEQ.start;
    return LOCSEQ.seq.substr(i,3).toUpperCase() === "ATG";
  }
  function baseAt(pos){ var i = pos - LOCSEQ.start; return (i>=0 && i<LOCSEQ.seq.length) ? LOCSEQ.seq.charAt(i).toUpperCase() : "N"; }

  function init(){
    if (META.palette){ for (var k in META.palette) if (PAL[k]!==undefined) PAL[k]=META.palette[k]; }
    TX_ORDER.forEach(function(t){ state[t]={open:false, sharedOpen:false}; });
    deriveIntrons(); deriveGPos(); buildLegendData();
    buildHeader(); buildLegend(); wireChrome();
    svg = d3s("stack"); hdr = d3s("hdr"); ideo = d3s("ideo");
    zt = d3.zoomIdentity;
    setupZoom();
    setupTouch();
    renderFeatureList();
    wireLegendToggle();
    // phone starts with side sheets CLOSED so the stack gets the full width (details-on-demand)
    setLeft(!isPhone()); openRight(false);
    rebuild();
    window.addEventListener("resize", onResizeDebounced);
  }
  // legend collapses behind its title on phone (saves vertical space); expanded elsewhere
  function wireLegendToggle(){
    var lg=document.querySelector(".legend"), tt=document.querySelector(".legend-title");
    if (!lg || !tt) return;
    if (isPhone()) lg.classList.add("collapsed");
    tt.addEventListener("click", function(){ lg.classList.toggle("collapsed"); });
  }

  // ---- derivations from meta ----
  function carriers(id){ return META.orfs[id].carriers; }
  function carrierCount(id){ return META.orfs[id].carriers.length; }
  function orfIdsOfTx(t){ return Object.keys(META.orfs).filter(function(id){ return carriers(id).indexOf(t)>=0; }); }
  function orfsOfClass(cls){ return Object.keys(META.orfs).filter(function(id){ return META.orfs[id].cls===cls; })
      .sort(function(a,b){ return (+a)-(+b); }); }
  function txCoding(t){ var c=META.transcripts[t].cds; return (c&&c.length)?"coding":"non-coding"; }
  function txSpan(t){ var e=META.transcripts[t].exons; return [e[0][0], e[e.length-1][1]]; }

  function deriveIntrons(){
    var ex = META.transcripts["201"].exons;
    CANON_INTRONS = [];
    for (var i=0;i<ex.length-1;i++) CANON_INTRONS.push([ex[i][1]+1, ex[i+1][0]-1]);
    TX_ORDER.forEach(function(t){
      RETAINED[t] = [];
      META.transcripts[t].exons.forEach(function(e){
        CANON_INTRONS.forEach(function(iv){ if (e[0]<iv[0] && e[1]>iv[1]) RETAINED[t].push(iv.slice()); });
      });
      var cds = META.transcripts[t].cds || [], cum = 0;
      CDS_PHASE[t] = cds.map(function(c){ var ph=(3-(cum%3))%3; cum += (c[1]-c[0]+1); return ph; });
    });
  }
  // spliced coding position -> genomic position, per ORF (from blocks; + strand ascending)
  function deriveGPos(){
    Object.keys(META.orfs).forEach(function(id){
      var g=[]; META.orfs[id].blocks.forEach(function(b){ for (var p=b[0]; p<=b[1]; p++) g.push(p); });
      GPOS[id]=g;
    });
  }
  // two-tier specificity (computed live from carriers): unique = 1 carrier, shared = >= 2
  function specGroups(t){
    var mem = orfIdsOfTx(t), uniq=[], shared=[];
    mem.forEach(function(id){ if (carrierCount(id)===1) uniq.push(id); else shared.push(id); });
    var byC = function(a,b){ return (carrierCount(a)-carrierCount(b)) || ((+a)-(+b)); };
    uniq.sort(byC); shared.sort(byC);
    return { uniq:uniq, shared:shared };
  }
  function exonSegments(t, exon){
    var cds = META.transcripts[t].cds || [], ri = RETAINED[t]||[];
    var cuts = {}; cuts[exon[0]]=1; cuts[exon[1]+1]=1;
    function addCuts(ivs){ ivs.forEach(function(iv){
      var a=Math.max(iv[0],exon[0]), b=Math.min(iv[1],exon[1]); if (a<=b){ cuts[a]=1; cuts[b+1]=1; } }); }
    addCuts(cds); addCuts(ri);
    var pts = Object.keys(cuts).map(Number).sort(function(a,b){return a-b;}), segs=[];
    for (var i=0;i<pts.length-1;i++){
      var a=pts[i], b=pts[i+1]-1; if (b<a) continue;
      var mid=(a+b)/2, type="utr";
      if (inAny(cds,mid)) type="cds"; else if (inAny(ri,mid)) type="ri";
      segs.push({a:a,b:b,type:type});
    }
    return segs;
  }
  function inAny(ivs,p){ for (var i=0;i<ivs.length;i++) if (p>=ivs[i][0] && p<=ivs[i][1]) return true; return false; }

  // ================= header / legend / chrome =================
  function buildHeader(){
    var nT = TX_ORDER.length, nO = Object.keys(META.orfs).length;
    $("tb-head").textContent = "MYCN locus";
    $("tb-sub").textContent = "— " + nT + " isoforms, " + nO + " encoded ORFs";
  }
  function buildLegendData(){
    LEGEND = [
      [PAL.nmyc,"N-Myc-frame CDS",false], [PAL.dmycn,"ΔMYCN (202) CDS",false],
      [PAL.ri_line,"206 retained intron",true], [PAL.utr,"UTR / non-coding",true],
      [PAL.pred_ext,"predicted N-term ext",false], [PAL.mycnot,"MYCNOT uORF",false],
      [PAL.musep,"MUSEP uORF",false], [PAL.term,"novel ORF",false]
    ];
  }
  function buildLegend(){   // T2.5: uniform swatches (no per-entry 'thin' class)
    $("legend-list").innerHTML = LEGEND.map(function(p){
      return '<li><i style="background:'+p[0]+'"></i>'+esc(p[1])+'</li>';
    }).join("");
  }
  function wireChrome(){
    var cx=$("cue-dismiss"); if (cx) cx.addEventListener("click", function(){ $("cue").hidden=true; });
    $("expand-all").addEventListener("click", function(){   // FIX E: fully expand every transcript (unique + shared)
      TX_ORDER.forEach(function(t){ state[t].open=true; state[t].sharedOpen=true; });   // T1.1
      rebuild(); renderFeatureList();
    });
    $("collapse-all").addEventListener("click", function(){
      TX_ORDER.forEach(function(t){ state[t].open=false; state[t].sharedOpen=false; });
      selectedTx=null;
      pickOrf(null);            // FIX 3: also clear any ORF selection + close Details (single-writer); this rebuilds + re-fits
      renderFeatureList();
    });
    $("download-svg").addEventListener("click", downloadSVG);
    $("download-png").addEventListener("click", downloadPNG);
    // control bar navigation
    $("locus-form").addEventListener("submit", function(ev){ ev.preventDefault(); gotoLocusInput($("locus-input").value); });
    $("zoom-in").addEventListener("click", function(){ zoomButton(1.6); });
    $("zoom-out").addEventListener("click", function(){ zoomButton(1/1.6); });
    $("reset-view").addEventListener("click", resetView);
    // lens + search + drawers
    $("lens-tx").addEventListener("click", function(){ setLens("tx"); });
    $("lens-orf").addEventListener("click", function(){ setLens("orf"); });
    $("feature-search").addEventListener("input", function(){ filterList(this.value); });
    // LEFT drawer pushes -> collapsing/opening changes viewer width -> relayout at the new width, keeping the view
    $("features-strip").addEventListener("click", function(){ var v=viewG.slice(); setLeft(true); rebuild(); applyRange(v[0],v[1],false); });
    $("features-collapse").addEventListener("click", function(){ var v=viewG.slice(); setLeft(false); rebuild(); applyRange(v[0],v[1],false); });
    $("details-strip").addEventListener("click", function(){ if (railHasContent) openRight(true); });
    // FIX 1 (desync): closing the panel DESELECTS via the single writer, so the variable and DOM never diverge
    $("details-close").addEventListener("click", function(){ pickOrf(null); });
  }

  // ================= drawers =================
  // LEFT (Features) pushes: collapsing/opening changes the viewer width -> rebuild() re-lays out.
  // On phone the drawers are full-screen sheets, so only ONE may be open at a time (opening one closes
  // the other). The guards can't recurse: setLeft only calls openRight when opening, and vice versa.
  function setLeft(open){ if(open && isPhone() && selectedOrf!=null) pickOrf(null);   // phone: opening Features deselects (keeps panel-visibility === selection, no desync)
    $("features").hidden=!open; $("features-strip").hidden=open; }
  // RIGHT (Details): ONE explicit source of truth. Always writes both DOM states to match `open`
  // (never a blind toggle), so every selection path lands in the same, correct state.
  function openRight(open){ rightOpen=!!open; if(rightOpen && isPhone()) setLeft(false); $("details").hidden=!rightOpen; $("details-strip").hidden=rightOpen; }

  // ================= zoom / navigation =================
  function plotL(){ return GUTTER; }
  function plotR(){ return W-PADR; }
  function setupZoom(){
    // Standard genome-browser / slippy-map convention (igv, UCSC, Google Maps):
    //   plain scroll wheel = ZOOM (cursor-centred, up=in/down=out)
    //   drag = PAN: horizontal drag pans the genomic X view; vertical drag scrolls the stack (Y)
    // The ruler (hdr) gets a brush (drag-select-to-zoom) — distinct target, no conflict.
    zoom = d3.zoom()
      .clickDistance(6)                    // moves < 6px still count as a click -> selection fires reliably
      .touchable(function(){ return false; })   // touch is handled by setupTouch() (pan+scroll+pinch), not d3
      .filter(zoomFilter)
      .on("start", function(ev){ dragPrevY = (ev.sourceEvent && ev.sourceEvent.clientY!=null) ? ev.sourceEvent.clientY : null; })
      .on("zoom", onZoom)
      // T1.2: after each gesture write the sanitized (y=0) transform back so internal __zoom.y
      // can never accumulate unbounded across vertical drags.
      .on("end", function(){ dragPrevY = null; try{ svg.property("__zoom", zt); }catch(e){} });
    svg.call(zoom);
  }
  // ---- touch gestures (the touch equivalents of the mouse scheme) ----
  //   ONE finger  = pan (horizontal genomic pan + vertical stack scroll), like a mouse drag
  //   TWO fingers = pinch-zoom the genomic X axis, centred on the pinch midpoint (like scroll-zoom)
  //   TAP (no move past threshold) = falls through to the native emulated click -> the SVG element's
  //     own click handler fires (expand transcript / select ORF). preventDefault only once a gesture
  //     starts, so a tap is never swallowed and a pan never fires a select.
  // d3-zoom's own touch is disabled (touchable=false) so these are the single touch code-path.
  var TCH = null, TAP_SLOP = 8;
  function svgScaleX(){ var r=svg.node().getBoundingClientRect(); return r.width ? (W/r.width) : 1; }
  function svgUserX(clientX){ var r=svg.node().getBoundingClientRect(); return (clientX-r.left)*(r.width?(W/r.width):1); }
  function setupTouch(){
    var node=svg.node();
    node.addEventListener("touchstart", tStart, {passive:false});
    node.addEventListener("touchmove",  tMove,  {passive:false});
    node.addEventListener("touchend",   tEnd,   {passive:false});
    node.addEventListener("touchcancel",tEnd,   {passive:false});
  }
  function tStart(ev){
    if (ev.touches.length===1){
      var t=ev.touches[0], sc=$("stack-scroll");
      TCH={ mode:"one", x0:t.clientX, y0:t.clientY, view0:viewG.slice(),
            st0:sc?sc.scrollTop:0, moved:false };
    } else if (ev.touches.length>=2){
      var a=ev.touches[0], b=ev.touches[1], mid=(a.clientX+b.clientX)/2;
      TCH={ mode:"pinch", dist0:Math.max(1,Math.hypot(a.clientX-b.clientX,a.clientY-b.clientY)),
            view0:viewG.slice(), anchorG:x.invert(svgUserX(mid)), moved:true };
      ev.preventDefault();   // two-finger gesture is unambiguously a pinch -> own it immediately
    }
  }
  function tMove(ev){
    if (!TCH) return;
    if (TCH.mode==="one" && ev.touches.length===1){
      var t=ev.touches[0], dx=t.clientX-TCH.x0, dy=t.clientY-TCH.y0;
      if (!TCH.moved && Math.hypot(dx,dy) < TAP_SLOP) return;   // still a tap candidate: let native click stand
      TCH.moved=true; ev.preventDefault();
      var wbp=TCH.view0[1]-TCH.view0[0];
      var bpPerUser=wbp/Math.max(1,(plotR()-plotL()));
      var a=TCH.view0[0] - dx*svgScaleX()*bpPerUser;         // horizontal genomic pan (absolute from start)
      applyRange(a, a+wbp, false);
      var sc=$("stack-scroll"); if(sc) sc.scrollTop=Math.max(0, TCH.st0 - dy);   // vertical stack scroll
    } else if (TCH.mode==="pinch" && ev.touches.length>=2){
      ev.preventDefault();
      var a2=ev.touches[0], b2=ev.touches[1], mid=(a2.clientX+b2.clientX)/2;
      var dist=Math.max(1,Math.hypot(a2.clientX-b2.clientX,a2.clientY-b2.clientY));
      var wbp=clamp((TCH.view0[1]-TCH.view0[0])/(dist/TCH.dist0), 20, DMAX-DMIN);
      var f=clamp((svgUserX(mid)-plotL())/Math.max(1,(plotR()-plotL())),0,1);
      var a=TCH.anchorG - f*wbp;                              // keep the anchored base under the pinch midpoint
      applyRange(a, a+wbp, false);
    }
  }
  function tEnd(ev){
    if (TCH && TCH.mode==="pinch" && ev.touches.length===1){
      // lifted to one finger mid-pinch: reseed a pan from the remaining finger so there's no jump
      var t=ev.touches[0], sc=$("stack-scroll");
      TCH={ mode:"one", x0:t.clientX, y0:t.clientY, view0:viewG.slice(), st0:sc?sc.scrollTop:0, moved:true };
      return;
    }
    if (ev.touches.length===0) TCH=null;   // a no-move tap emulates a click -> element handler selects/expands
  }

  function zoomFilter(ev){
    if (ev.type === "wheel") return true;      // plain wheel zooms (cursor-centred)
    if (ev.type === "dblclick") return false;
    return !ev.button;                         // left-drag pans (x) + scrolls the stack (y)
  }
  // drag down (clientY increases -> dy>0) reveals rows above -> scrollTop decreases
  function scrollStackBy(dy){ var sc=$("stack-scroll"); if(!sc) return 0; sc.scrollTop = Math.max(0, sc.scrollTop - dy); return sc.scrollTop; }
  function refreshZoomExtent(){
    var maxk = (DMAX-DMIN) / 40;            // deepest ~40 bp visible
    // explicit extent (don't probe the DOM node's viewBox) + constrain to the plot area
    zoom.extent([[plotL(),0],[plotR(),Math.max(1,H)]])
        .scaleExtent([1, Math.max(2, maxk)])
        .translateExtent([[plotL(),-Infinity],[plotR(),Infinity]]);
    svg.call(zoom.transform, zt);
  }
  // Zoom stays constrained to the genomic X axis (rows are never scaled). A vertical DRAG
  // component is translated into native scroll of the stack container (move up/down the stack).
  function onZoom(ev){
    zt = d3.zoomIdentity.translate(ev.transform.x,0).scale(ev.transform.k);
    var se = ev.sourceEvent;
    if (se && dragPrevY!=null && se.clientY!=null && (se.type==="mousemove"||se.type==="pointermove")){
      var dy = se.clientY - dragPrevY; dragPrevY = se.clientY;
      if (dy) scrollStackBy(dy);   // vertical drag -> scroll the stack (rows move up/down)
    }
    scheduleDraw();   // T1.5: coalesce many zoom/mousemove events into one render per animation frame
  }
  // T1.5: rAF-throttled redraw. Headless tests (NOANIM) render synchronously for determinism.
  var drawPending=false, RAF=(typeof requestAnimationFrame==="function")?requestAnimationFrame:function(f){return setTimeout(f,16);};
  function scheduleDraw(){
    if (NOANIM){ draw(); return; }
    if (drawPending) return;
    drawPending=true; RAF(function(){ drawPending=false; draw(); });
  }
  function zoomButton(factor){
    var c = (plotL()+plotR())/2;
    svg.transition().duration(180).call(zoom.scaleBy, factor, [c,0]);
  }
  // navigate so [g0,g1] fills the plot (feature-padded)
  function flyTo(g0, g1, animate){
    var pad = Math.max(30, (g1-g0)*0.12);
    applyRange(clamp(g0-pad,DMIN,DMAX), clamp(g1+pad,DMIN,DMAX), animate);
  }
  // navigate to EXACTLY [g0,g1] (tight; used by ruler-brush + locus box so narrow spans reach base level)
  function zoomToRange(g0, g1, animate){
    if (g1<g0){ var s=g0; g0=g1; g1=s; }
    applyRange(clamp(g0,DMIN,DMAX), clamp(g1,DMIN,DMAX), animate);
  }
  function applyRange(a, b, animate){
    if (b-a < 20){ var m=(a+b)/2; a=clamp(m-10,DMIN,DMAX); b=clamp(m+10,DMIN,DMAX); }  // floor ~20 bp
    var k = clamp((plotR()-plotL())/(x0(b)-x0(a)), 1, zoom.scaleExtent()[1]);
    var tx = plotL() - k*x0(a);
    var t = d3.zoomIdentity.translate(tx,0).scale(k);
    if (animate && !NOANIM) svg.transition().duration(360).call(zoom.transform, t);
    else svg.call(zoom.transform, t);
  }
  var NOANIM = false;   // set true only by headless tests so transitions apply synchronously
  function resetView(){
    selectedTx=null;
    pickOrf(null);              // single-writer clear of the ORF selection (closes Details, clears highlight/dim)
    flyTo(DMIN, DMAX, true);    // then reset the view to the full locus
    setNow("MYCN locus", null);
    renderFeatureList();
  }
  function gotoLocusInput(str){
    var m = String(str).replace(/,/g,"").match(/(?:chr2[:\s])?\s*(\d+)\s*[-–:]\s*(\d+)/i);
    if (!m){                                    // T2.3: give feedback on unparseable input
      var el=$("locus-input"); if(el){ el.classList.add("bad"); setTimeout(function(){ el.classList.remove("bad"); },900); }
      return;
    }
    var a=+m[1], b=+m[2]; if (a>b){ var t=a;a=b;b=t; }
    zoomToRange(clamp(a,DMIN,DMAX), clamp(b,DMIN,DMAX), true);
  }

  // ================= layout =================
  // device class from the VIEWPORT width (what "phone/tablet/desktop" means), independent of drawer state
  function vpW(){ return (typeof window!=="undefined" && window.innerWidth) ? window.innerWidth : 1024; }
  function isPhone(){ return vpW() < 640; }
  function isTablet(){ var w=vpW(); return w>=640 && w<1024; }
  // GUTTER (label lane) + PADR shrink on small screens so the plot isn't crowded; W tracks the container
  function computeGutter(){
    if (isPhone()){ GUTTER=92; PADR=8; }
    else if (isTablet()){ GUTTER=150; PADR=12; }
    else { GUTTER=200; PADR=18; }
  }
  function computeW(){
    computeGutter();
    var box=$("stack-scroll"), cw=(box?box.clientWidth:900);
    var floor = isPhone()?300:(isTablet()?480:720);   // desktop/jsdom (cw=0) -> 720, as before
    W = Math.max(floor, cw-2);
  }
  var SHORT_LBL_MAX = 130;   // GUTTER below this -> shorten transcript labels ("201", not "MYCN-201")
  function shortLabels(){ return GUTTER < SHORT_LBL_MAX; }
  function layout(){
    var rows=[], y=8;
    TX_ORDER.forEach(function(t){
      rows.push({kind:"tx", t:t, y:y}); y+=TXH;
      if (state[t].open){
        var g=specGroups(t);
        y=pushGroup(rows,"Unique to MYCN-"+t,g.uniq,true,t,y);
        if (g.shared.length){
          rows.push({kind:"shared",t:t,y:y,n:g.shared.length}); y+=SHH;
          if (state[t].sharedOpen) g.shared.forEach(function(id){ rows.push({kind:"orf",t:t,id:id,uniq:false,y:y}); y+=ORFH; });
        }
        y+=EXPAD;
      }
    });
    H = y+12;
    return rows;
  }
  function pushGroup(rows,label,ids,uniq,t,y){
    if (!ids.length) return y;
    rows.push({kind:"grp",label:label,uniq:uniq,y:y}); y+=GRPH;
    ids.forEach(function(id){ rows.push({kind:"orf",t:t,id:id,uniq:uniq,y:y}); y+=ORFH; });
    return y;
  }

  // ================= render orchestration =================
  function rebuild(){
    computeW();
    x0 = d3.scaleLinear().domain([DMIN,DMAX]).range([plotL(),plotR()]);
    ROWS = layout();
    svg.attr("width",W).attr("height",H).attr("viewBox","0 0 "+W+" "+H).attr("font-family",FONT);
    hdr.attr("width",W).attr("height",HDR_H).attr("viewBox","0 0 "+W+" "+HDR_H).attr("font-family",FONT);
    ideo.attr("width",W).attr("height",IDEO_H).attr("viewBox","0 0 "+W+" "+IDEO_H).attr("font-family",FONT);
    refreshZoomExtent();
    drawIdeogramBar();
    draw();
  }
  function onResize(){ var v=viewG.slice(); rebuild(); applyRange(v[0],v[1],false); }
  // debounced resize: recompute GUTTER + W, relayout, and redraw at the new width, keeping the current view
  var resizeTimer=null;
  function onResizeDebounced(){ if(resizeTimer) clearTimeout(resizeTimer); resizeTimer=setTimeout(onResize,120); }

  function draw(){
    x = zt.rescaleX(x0);
    viewG = [clamp(x.invert(plotL()),DMIN,DMAX), clamp(x.invert(plotR()),DMIN,DMAX)];
    renderHeader();
    renderStack();
    updateIdeoMarker();
    updateNow();
    $("locus-input").value = "chr2:"+fmt(viewG[0])+"-"+fmt(viewG[1]);
  }

  // ---- chevrons (rightward, + strand) ----
  function chevrons(g, xa, xb, mid, color){
    for (var cx=xa+26; cx<xb-5; cx+=52){
      g.append("path").attr("d","M"+(cx-3)+" "+(mid-4)+" L"+(cx+3)+" "+mid+" L"+(cx-3)+" "+(mid+4))
        .attr("fill","none").attr("stroke",color).attr("stroke-width",1.1);
    }
  }
  function pxPerBase(){ return (plotR()-plotL())/(viewG[1]-viewG[0]); }

  // ================= chromosome ideogram (real cytobands; FIXED chromosome-scale) =================
  // The ideogram is a whole-chr2 overview at its own chromosome scale (ideoX), entirely
  // independent of the locus x-scale. Only the locus marker moves; the banded bar is static.
  function ideoGeom(){
    var chrLen = (CYTO && CYTO.length) || CHR2LEN;
    var L = GUTTER, R = W-PADR, midY = 26, barH = 14;
    var ix = function(p){ return L + (R-L)*(clamp(p,0,chrLen)/chrLen); };
    return { chrLen:chrLen, L:L, R:R, midY:midY, barH:barH, ix:ix };
  }
  function acenSpan(){
    var a=(CYTO.bands||[]).filter(function(b){ return b.stain==="acen"; });
    if (!a.length) return null;
    return { start:a[0].start, mid:a[0].end, end:a[a.length-1].end };
  }
  function bandAt(pos){
    var bs=(CYTO&&CYTO.bands)||[];
    for (var i=0;i<bs.length;i++) if (pos>=bs[i].start && pos<bs[i].end) return bs[i];
    return null;
  }
  function drawIdeogramBar(){
    ideo.selectAll("*").remove();
    if (!CYTO || !CYTO.bands){ return; }
    var G=ideoGeom(), ix=G.ix, midY=G.midY, barH=G.barH, top=midY-barH/2, bot=midY+barH/2;
    var acen=acenSpan();
    ideo.append("text").attr("x",12).attr("y",midY+4).attr("font-size",11).attr("font-weight",600)
      .attr("fill","#5f6a66").text("chr2");
    // rounded-chromosome clip so band ends look like a chromosome
    var defs=ideo.append("defs");
    defs.append("clipPath").attr("id","ideo-clip").append("rect")
      .attr("x",ix(0)).attr("y",top).attr("width",ix(G.chrLen)-ix(0)).attr("height",barH).attr("rx",barH/2).attr("ry",barH/2);
    var gb=ideo.append("g").attr("clip-path","url(#ideo-clip)");
    // non-acen bands as scaled rects; acen handled as a pinched constriction
    CYTO.bands.forEach(function(b){
      if (b.stain==="acen") return;
      var xa=ix(b.start), xb=ix(b.end);
      gb.append("rect").attr("x",xa).attr("y",top).attr("width",Math.max(0.4,xb-xa)).attr("height",barH)
        .attr("fill", STAIN[b.stain]||"#d9dee2")
        .append("title").text("chr2 "+b.name+" ("+b.stain+")");
    });
    // centromere: hourglass/bowtie meeting at the constriction (real acen coords)
    if (acen){
      var xs=ix(acen.start), xm=ix(acen.mid), xe=ix(acen.end);
      gb.append("path").attr("d","M"+xs+" "+top+" L"+xm+" "+midY+" L"+xe+" "+top+" Z").attr("fill",STAIN.acen);
      gb.append("path").attr("d","M"+xs+" "+bot+" L"+xm+" "+midY+" L"+xe+" "+bot+" Z").attr("fill",STAIN.acen);
    }
    // chromosome outline
    ideo.append("rect").attr("x",ix(0)).attr("y",top).attr("width",ix(G.chrLen)-ix(0)).attr("height",barH)
      .attr("rx",barH/2).attr("ry",barH/2).attr("fill","none").attr("stroke","#9aa4ab").attr("stroke-width",1);
    // locus marker (updated by updateIdeoMarker)
    ideo.append("path").attr("id","ideo-mark-tip").attr("fill","#c0392b");
    ideo.append("rect").attr("id","ideo-mark").attr("y",top-3).attr("height",barH+6).attr("rx",2)
      .attr("fill","none").attr("stroke","#c0392b").attr("stroke-width",1.7);
    ideo.append("text").attr("id","ideo-band").attr("y",bot+13).attr("font-size",9.5).attr("font-weight",600)
      .attr("fill","#a5322c").attr("text-anchor","middle").text("");
    updateIdeoMarker();
  }
  function updateIdeoMarker(){
    if (!CYTO || !CYTO.bands) return;
    var G=ideoGeom(), ix=G.ix, midY=G.midY, barH=G.barH, top=midY-barH/2;
    var mid=(viewG[0]+viewG[1])/2;
    var xa=ix(viewG[0]), xb=ix(viewG[1]), cx=ix(mid);
    var w=Math.max(6, xb-xa);                 // locus is sub-pixel on the whole chromosome -> min width
    var m=ideo.select("#ideo-mark"); if(m.empty()) return;
    m.attr("x",cx-w/2).attr("width",w);
    ideo.select("#ideo-mark-tip").attr("d","M"+(cx-4)+" "+(top-6)+" L"+(cx+4)+" "+(top-6)+" L"+cx+" "+(top-1)+" Z");
    var band=bandAt(mid), name=band?band.name:"";
    // clamp the (middle-anchored) label so the whole word stays inside the SVG width — never
    // runs off the left/right edge when the marker sits near a chromosome end.
    var halfW = name.length*0.62*9.5/2 + 3;
    var lx = clamp(cx, halfW+2, W-halfW-2);
    ideo.select("#ideo-band").attr("x",lx).text(name);
  }

  // ================= genomic header: ruler + reference strip + codon track =================
  function renderHeader(){
    hdr.selectAll("*").remove();
    hdr.append("defs").append("clipPath").attr("id","hdr-clip").append("rect")
      .attr("x",plotL()).attr("y",0).attr("width",Math.max(1,plotR()-plotL())).attr("height",HDR_H);
    var rulerY=14, stripY=24, stripH=13, codonY=stripY+stripH+2;
    // opaque gutter cover so ruler/ref/base letters never show under the "ref" label lane
    hdr.append("rect").attr("x",0).attr("y",0).attr("width",plotL()-1).attr("height",HDR_H).attr("fill","#ffffff");
    // ruler (clipped to the plot)
    var g=hdr.append("g").attr("clip-path","url(#hdr-clip)");
    g.append("line").attr("x1",plotL()).attr("x2",plotR()).attr("y1",rulerY).attr("y2",rulerY).attr("stroke","#d7dde1");
    // FIX 3: responsive tick density — target count from plot width, then greedily drop any label that
    // would come within 8px of the previous one (measured), so labels are never a mashed-together smear.
    var plotW=plotR()-plotL();
    var target=Math.max(2, Math.min(8, Math.round(plotW/90)));   // ~90px per label incl. gap
    var lastRight=-1e9;
    x.ticks(target).forEach(function(v){
      var xp=x(v); if (xp<plotL()-1||xp>plotR()+1) return;
      g.append("line").attr("x1",xp).attr("x2",xp).attr("y1",rulerY-4).attr("y2",rulerY).attr("stroke","#c4ccd1");
      var txt=g.append("text").attr("class","ruler-lbl").attr("x",xp).attr("y",rulerY-6).attr("text-anchor","middle")
        .attr("font-size",9.5).attr("font-family",MONO).attr("fill","#5f6a66")
        .text((v/1000).toFixed(2)+" kb");
      var w=(txt.node().getComputedTextLength?txt.node().getComputedTextLength():0);
      if (w && (xp - w/2) < lastRight + 8) txt.remove();          // would overlap the previous label -> drop it
      else lastRight = xp + (w||0)/2;
    });
    hdr.append("text").attr("x",12).attr("y",stripY+stripH-3).attr("font-size",9.5).attr("font-weight",600)
      .attr("fill","#5f6a66").text("ref");
    renderReferenceStrip(hdr.append("g").attr("clip-path","url(#hdr-clip)"), stripY, stripH);
    renderCodonTrack(hdr.append("g").attr("clip-path","url(#hdr-clip)"), codonY);
    installRulerBrush(codonY+13);
  }
  // FIX 2: drag across the ruler/ref strip to select a span, release to zoom to it (igv behaviour).
  // Distinct target from the stack body: brush on hdr, pan-zoom on svg -> no gesture conflict.
  function installRulerBrush(yBottom){
    rulerBrush = d3.brushX()
      .extent([[plotL(),1],[plotR(),Math.max(20,yBottom)]])
      .on("end", onRulerBrushEnd);
    var gB = hdr.append("g").attr("class","ruler-brush").call(rulerBrush);
    gB.select(".overlay").attr("cursor","ew-resize");
  }
  function onRulerBrushEnd(ev){
    if (!ev.selection) return;                          // plain click (no drag) -> ignore
    var px0=ev.selection[0], px1=ev.selection[1];
    hdr.select(".ruler-brush").call(rulerBrush.move, null);   // clear the selection rectangle
    if (px1-px0 < 4) return;                            // guard against zero/near-zero drags
    zoomToRange(x.invert(px0), x.invert(px1), true);    // reuse the existing x.invert + zoom path
  }
  function renderReferenceStrip(g, y, h){
    var ppb=pxPerBase();
    var a=Math.max(DMIN, Math.floor(viewG[0])), b=Math.min(DMAX, Math.ceil(viewG[1]));
    if (ppb >= LETTER_PXPB){
      for (var p=a; p<=b; p++){
        var base=baseAt(p), xa=x(p), xb=x(p+1), cx=(xa+xb)/2;
        g.append("rect").attr("x",xa).attr("y",y).attr("width",Math.max(1,xb-xa)).attr("height",h)
          .attr("fill",BASECOL[base]||BASECOL.N).attr("opacity",.20);
        g.append("text").attr("class","base-let").attr("x",cx).attr("y",y+h-2.5).attr("text-anchor","middle")
          .attr("font-size",Math.min(13,(xb-xa)*0.9)).attr("font-family",MONO).attr("font-weight",600)
          .attr("fill",BASECOL[base]||BASECOL.N).text(base);
      }
    } else {
      var step = ppb>=CELL_PXPB ? 1 : Math.max(1, Math.ceil((b-a)/800));
      for (var q=a; q<=b; q+=step){
        var bb=baseAt(q), x1=x(q), x2=x(q+step);
        g.append("rect").attr("x",x1).attr("y",y).attr("width",Math.max(0.6,x2-x1)).attr("height",h)
          .attr("fill",BASECOL[bb]||BASECOL.N).attr("opacity", ppb>=CELL_PXPB?0.9:0.75);
      }
    }
  }
  // codon/aa track for the SELECTED ORF (spliced->genomic map handles intron-spanning codons correctly)
  function renderCodonTrack(g, y){
    if (pxPerBase() < LETTER_PXPB || !selectedOrf) return;
    var id=selectedOrf, o=META.orfs[id], gp=GPOS[id], aa=o.aa||"", nt=o.nt||"";
    var ncod=Math.floor(gp.length/3);
    for (var k=0;k<ncod;k++){
      var p0=gp[3*k], p1=gp[3*k+1], p2=gp[3*k+2];
      if (Math.max(p0,p1,p2) < viewG[0]-2 || Math.min(p0,p1,p2) > viewG[1]+2) continue;
      var res = k<aa.length ? aa.charAt(k) : "*";
      var isStart=(k===0), isStop=(k>=aa.length);
      var fill = isStart?"#1f6b45":(isStop?"#a5322c":"#5b6b82");
      var bg   = isStart?"#d9efe1":(isStop?"#f3dcda":"#eef2f7");
      // T1.6: anchor the aa over the codon's bases that are IN VIEW (exonic/on-screen), not the
      // p0..p2 midpoint — for a junction-spanning codon the midpoint falls in the intron gap.
      var inView=[p0,p1,p2].filter(function(pp){ return pp>=viewG[0]-2 && pp<=viewG[1]+2; }).map(function(pp){ return x(pp+0.5); });
      var cx = inView.length ? inView.reduce(function(s,v){return s+v;},0)/inView.length : x(p1+0.5);
      // underline spanning the 3 bases (crosses the intron honestly if spliced)
      g.append("line").attr("x1",x(Math.min(p0,p2)+0)).attr("x2",x(Math.max(p0,p2)+1))
        .attr("y1",y+11).attr("y2",y+11).attr("stroke",fill).attr("stroke-width",1).attr("opacity",.5);
      g.append("rect").attr("x",cx-5).attr("y",y).attr("width",10).attr("height",12).attr("rx",2).attr("fill",bg);
      g.append("text").attr("x",cx).attr("y",y+9.5).attr("text-anchor","middle").attr("font-size",9.5)
        .attr("font-family",MONO).attr("font-weight",600).attr("fill",fill).text(res);
    }
    hdr.append("text").attr("x",12).attr("y",y+10).attr("font-size",9).attr("font-weight",600).attr("fill","#5f6a66")
      .text("ORF"+id);   // track label in the gutter lane (unclipped)
  }

  // ================= stack rows =================
  // FIX A root cause: genomic STRUCTURE is drawn inside clip-path #plot-clip (x GUTTER..W-PADR),
  // so exon/intron/CDS/ORF blocks + chevrons never paint into the label gutter on zoom/pan.
  // Labels live UNCLIPPED in the gutter lane, over an opaque gutter background.
  function isCarrierOf(orfId,t){ return !!orfId && META.orfs[orfId].carriers.indexOf(t)>=0; }
  function rowOpacityFor(t){ return (selectedOrf && !isCarrierOf(selectedOrf,t)) ? 0.35 : 1; }   // FIX D dim
  function clipped(g){ return g.append("g").attr("clip-path","url(#plot-clip)"); }

  function renderStack(){
    svg.selectAll("*").remove();
    svg.append("defs").append("clipPath").attr("id","plot-clip").append("rect")
      .attr("x",plotL()).attr("y",0).attr("width",Math.max(1,plotR()-plotL())).attr("height",H);
    gMain = svg.append("g");
    // FIX 1: empty-space deselect catcher at the BOTTOM — a tap on empty plot/gutter clears the ORF
    // selection. Decorative layers above are pointer-events:none so the click reaches this catcher;
    // feature hit-rects (drawn later, on top) still capture their own clicks. A drag (>clickDistance)
    // is a pan (d3-zoom on the svg), never a deselect.
    gMain.append("rect").attr("class","deselect-catcher").attr("x",0).attr("y",0).attr("width",W).attr("height",H)
      .attr("fill","transparent").on("click",function(){ pickOrf(null); });
    // opaque gutter background (belt-and-braces: nothing genomic ever shows under the labels)
    gMain.append("rect").attr("class","gutter-bg").attr("x",0).attr("y",0).attr("width",plotL()-1).attr("height",H)
      .attr("fill","#ffffff").style("pointer-events","none");
    // gridlines + selection band -> clipped to the plot; decorative -> pointer-events:none
    var gGrid=clipped(gMain).style("pointer-events","none");
    x.ticks(8).forEach(function(v){ var xp=x(v); if (xp<plotL()||xp>plotR()) return;
      gGrid.append("line").attr("x1",xp).attr("x2",xp).attr("y1",0).attr("y2",H).attr("stroke","#eef1f3"); });
    if (selectedOrf){ var sp=META.orfs[selectedOrf].span;
      gGrid.append("rect").attr("x",x(sp[0])).attr("y",0).attr("width",Math.max(1.5,x(sp[1]+1)-x(sp[0])))
        .attr("height",H).attr("fill",PAL.nmyc).attr("opacity",0.06); }
    // gutter divider (label lane, unclipped, decorative)
    gMain.append("line").attr("x1",plotL()-6).attr("x2",plotL()-6).attr("y1",0).attr("y2",H).attr("stroke","#e3e7e9").style("pointer-events","none");
    ROWS.forEach(function(r){
      if (r.kind==="tx") drawTx(r);
      else if (r.kind==="grp") drawGroup(r);
      else if (r.kind==="orf") drawOrf(r);
      else if (r.kind==="shared") drawShared(r);
    });
  }
  // T2.1: the ΔMYCN transcript is the one carrying ORF2 (ΔMYCN) — derived, not hardcoded "202"
  function dmycnTx(){ return (META.orfs["2"] && META.orfs["2"].carriers[0]) || "202"; }
  function segFill(t,type){ if(type==="cds") return (t===dmycnTx())?PAL.dmycn:PAL.nmyc; if(type==="ri") return PAL.ri_line; return PAL.utr; }

  // FIX D: emphasise the selected ORF's blocks on a carrier transcript's row (halo + solid + accent outline)
  function drawOrfHighlight(gClip, orfId, mid){
    var o=META.orfs[orfId];
    var b0=x(o.blocks[0][0]), b1=x(o.blocks[o.blocks.length-1][1]+1);
    gClip.append("line").attr("x1",b0).attr("x2",b1).attr("y1",mid).attr("y2",mid).attr("stroke",PAL.nmyc).attr("stroke-width",1).attr("opacity",.55).style("pointer-events","none");
    o.blocks.forEach(function(b){
      var xa=x(b[0]), xb=x(b[1]+1), w=Math.max(2,xb-xa);
      gClip.append("rect").attr("class","orf-hl-halo").attr("x",xa-1.5).attr("y",mid-8).attr("width",w+3).attr("height",16).attr("rx",3).attr("fill",o.color).attr("opacity",.30);
      // FIX 1(b): the highlighted ORF block on a carrier row is the OBVIOUS toggle target, so give it the ORF's
      // own click -> pickOrf(orfId), which CLOSES on a second click (orfId===selectedOrf here). Previously it was
      // a filled, handler-less rect painted OVER the exon: it swallowed the click, so the toggle never fired and
      // only a white-space click closed the panel. (.orf-hl-halo/.orf-sel/the connector line are pointer-events:
      // none decoration — see styles.css:157 — so only this .orf-hl block is clickable.)
      gClip.append("rect").attr("class","orf-hl").attr("x",xa).attr("y",mid-6).attr("width",w).attr("height",12).attr("rx",2.5)
        .attr("fill",o.color).attr("stroke",PAL.nmyc).attr("stroke-width",1.6)
        .on("click",function(ev){ ev.stopPropagation(); pickOrf(orfId); })
        .on("mouseenter",function(ev){ tip(ev,orfTip(orfId)); }).on("mousemove",moveTip).on("mouseleave",hideTip);
    });
  }
  // bring the first carrier transcript row into view (vertical scroll only; do not fly-zoom away)
  function scrollCarrierIntoView(orfId){
    var sc=$("stack-scroll"); if(!sc||!orfId) return;
    var cs=META.orfs[orfId].carriers;
    var row=null; for (var i=0;i<ROWS.length;i++){ if (ROWS[i].kind==="tx" && cs.indexOf(ROWS[i].t)>=0){ row=ROWS[i]; break; } }
    if (!row) return;
    var top=row.y, bottom=row.y+TXH, vt=sc.scrollTop, vh=sc.clientHeight||400;
    if (top < vt) sc.scrollTop = Math.max(0, top-8);
    else if (bottom > vt+vh) sc.scrollTop = bottom - vh + 8;
  }

  function drawTx(r){
    var t=r.t, tx=META.transcripts[t], mid=r.y+TXH/2, op=rowOpacityFor(t);
    var g=gMain.append("g").attr("class","tx-row").attr("data-tx",t).attr("opacity",op);   // T2.4: dropped vestigial .hot (no CSS; selection shown via label fill)
    var gc=clipped(g);   // genomic structure (clipped to the plot)
    // ---- labels (unclipped, gutter lane) ----
    g.append("rect").attr("class","gutter-hit").attr("x",0).attr("y",r.y).attr("width",GUTTER-8).attr("height",TXH)
      .attr("fill","transparent").on("click",function(){ pickTx(t); })
      .on("mouseenter",function(ev){ tip(ev,"<b>MYCN-"+t+"</b> · "+orfIdsOfTx(t).length+" ORFs<br><span class='dim'>click to "+(state[t].open?"collapse":"expand")+"</span>"); })
      .on("mousemove",moveTip).on("mouseleave",hideTip);
    // on a narrow (phone) gutter, shorten "MYCN-201" -> "201" and drop "coding ·" from the meta line;
    // the full name + ORF count stay reachable via the row tooltip (hover) and tap-to-expand.
    var sh=shortLabels(), nOrf=orfIdsOfTx(t).length;
    g.append("text").attr("class","tx-label").attr("x",14).attr("y",mid-3).attr("font-size",13).attr("font-weight",600)
      .attr("fill",selectedTx===t?PAL.nmyc:"#222a2e").text(sh ? t : "MYCN-"+t);
    var coding=txCoding(t);
    g.append("text").attr("class","tx-meta").attr("x",14).attr("y",mid+12).attr("font-size",10.5)
      .attr("fill",coding==="non-coding"?"#a5322c":"#5f6a66").text(sh ? (nOrf+" ORFs") : (coding+" · "+nOrf+" ORFs"));
    // right-edge expand chevron (row disclosure indicator; the whole gutter row is the click target via
    // gutter-hit -> pickTx). FIX 5: the per-row "zoom to transcript" magnifier was removed — zooming to one
    // of 8 near-identical ~9 kb isoforms barely differs and clips the others out of frame; the locus box,
    // +/-, scroll/pinch, ruler-brush and Reset cover every real navigation need.
    var zx=GUTTER-16;
    g.append("text").attr("class","tx-chev").attr("x",zx).attr("y",mid+4).attr("font-size",12)
      .attr("fill","#5f6a66").attr("text-anchor","middle").text(state[t].open?"▾":"▸");
    // ---- structure (clipped) ----
    var x0e=x(tx.exons[0][0]), x1e=x(tx.exons[tx.exons.length-1][1]+1);
    gc.append("line").attr("class","backbone").attr("x1",x0e).attr("x2",x1e).attr("y1",mid).attr("y2",mid)
      .attr("stroke",BACKBONE).attr("stroke-width",1.2);
    chevrons(gc,x0e,x1e,mid,BACKBONE);
    tx.exons.forEach(function(exon,ei){
      exonSegments(t,exon).forEach(function(s){
        var h=s.type==="cds"?THICK:(s.type==="ri"?RIH:THIN), xa=x(s.a), xb=x(s.b+1), w=Math.max(1.2,xb-xa);
        gc.append("rect").attr("class","exon").attr("x",xa).attr("y",mid-h/2).attr("width",w).attr("height",h).attr("rx",1.5)
          .attr("fill",segFill(t,s.type))
          .on("click",function(ev){ ev.stopPropagation(); showExon(t,exon,ei); })
          .on("mouseenter",function(ev){ tip(ev,exonTip(t,exon,ei,s)); }).on("mousemove",moveTip).on("mouseleave",hideTip);
      });
    });
    // FIX D: light up the selected ORF where it sits on THIS carrier transcript
    if (isCarrierOf(selectedOrf,t)) drawOrfHighlight(gc, selectedOrf, mid);
  }
  function drawGroup(r){
    var t=gMain.append("text").attr("class","grp-lbl").attr("x",20).attr("y",r.y+14).attr("font-size",10)
      .attr("font-weight",600).attr("letter-spacing",".5px").attr("fill",r.uniq?PAL.nmyc:"#5f6a66")
      .text(r.label.toUpperCase());
    if (r.uniq) t.append("title").text("within the 8 transcripts shown");   // FIX 2: scope, on hover
  }
  function drawShared(r){
    var t=r.t, g=gMain.append("g");
    g.append("rect").attr("class","shared-toggle-hit").attr("x",20).attr("y",r.y+3).attr("width",GUTTER-40).attr("height",SHH-8)
      .attr("rx",5).attr("fill","#f2f4f6").attr("stroke","#d5dbe2").attr("stroke-dasharray","3 3")
      .on("click",function(){ state[t].sharedOpen=!state[t].sharedOpen; rebuild(); });
    g.append("text").attr("class","shared-toggle-lbl").attr("x",30).attr("y",r.y+SHH/2+2).attr("font-size",11)
      .attr("font-weight",600).attr("fill",PAL.nmyc).text((state[t].sharedOpen?"▾ hide ":"▸ show ")+r.n+" shared ORFs");
  }
  function drawOrf(r){
    var id=r.id, o=META.orfs[id], mid=r.y+ORFH/2;
    var g=gMain.append("g").attr("class","orf-g").attr("data-id",id);
    var gc=clipped(g);
    // ---- label columns (gutter): swatch | ORF{id}{name} | aa·carrier (right-aligned) ----
    g.append("rect").attr("class","orf-hit").attr("x",22).attr("y",r.y).attr("width",GUTTER-30).attr("height",ORFH).attr("fill","transparent")
      .style("cursor","pointer").on("click",function(){ pickOrf(id); })
      .on("mouseenter",function(ev){ tip(ev,orfTip(id)); }).on("mousemove",moveTip).on("mouseleave",hideTip);
    g.append("rect").attr("class","orf-sw").attr("x",24).attr("y",mid-5).attr("width",10).attr("height",10).attr("rx",2).attr("fill",o.color);
    var nm=NAME[id]?(" "+NAME[id]):"", TAGX=GUTTER-8, LBLX=42;
    // narrow (phone) gutter: tag shows only carrier count so the id+name column keeps room before truncation
    var tagStr=shortLabels() ? (carrierCount(id)+"c") : (o.aa_len+"aa · "+carrierCount(id)+"c");
    var tag=g.append("text").attr("class","orf-tag").attr("x",TAGX).attr("y",mid+3.5).attr("text-anchor","end")
      .attr("font-size",10).attr("fill","#5f6a66").text(tagStr);
    var lbl=g.append("text").attr("class","orf-lbl").attr("x",LBLX).attr("y",mid+3.5).attr("font-size",11)
      .attr("fill","#222a2e").on("click",function(){ pickOrf(id); }).text("ORF"+id+nm);
    // measured truncation (browser): keep >= 8px between the id+name and the tag column
    try{
      var tn=tag.node(), ln=lbl.node();
      var tagW=(tn.getComputedTextLength?tn.getComputedTextLength():0);
      var avail=(TAGX - tagW - 8) - LBLX, full="ORF"+id+nm, s=full;
      while (ln.getComputedTextLength && ln.getComputedTextLength()>avail && s.length>5){ s=s.slice(0,-1); ln.textContent=s+"…"; }
    }catch(e){}
    // ---- structure (clipped) ----
    var b0=x(o.blocks[0][0]), b1=x(o.blocks[o.blocks.length-1][1]+1);
    gc.append("line").attr("x1",b0).attr("x2",b1).attr("y1",mid).attr("y2",mid).attr("stroke","#c2ccc9").attr("stroke-width",1);
    chevrons(gc,b0,b1,mid,"#c2ccc9");
    o.blocks.forEach(function(b){
      var xa=x(b[0]), xb=x(b[1]+1), w=Math.max(1.5,xb-xa);
      gc.append("rect").attr("class","orf-bar").attr("x",xa).attr("y",mid-ORF_THICK/2).attr("width",w).attr("height",ORF_THICK)
        .attr("rx",2).attr("fill",o.color).on("click",function(){ pickOrf(id); })
        .on("mouseenter",function(ev){ tip(ev,orfTip(id)); }).on("mousemove",moveTip).on("mouseleave",hideTip);
    });
    if (selectedOrf===id){
      gc.append("rect").attr("class","orf-sel").attr("x",b0-3).attr("y",mid-ORF_THICK/2-3).attr("width",(b1-b0)+6).attr("height",ORF_THICK+6)
        .attr("rx",4).attr("fill","none").attr("stroke",PAL.nmyc).attr("stroke-width",1.6);
    }
  }

  // ================= now-showing =================
  function setNow(name, g){ selNowName=name; drawNow(name, g||viewG); }
  var selNowName="MYCN locus";
  function updateNow(){ drawNow(selNowName, viewG); }
  function drawNow(name, g){
    $("ns-name").textContent = name;
    $("ns-coord").textContent = "chr2:"+fmt(g[0])+"-"+fmt(g[1]);
  }

  // ================= interactions =================
  function paddedSpan(sp){ var pad=Math.max(30,(sp[1]-sp[0])*0.12); return [clamp(sp[0]-pad,DMIN,DMAX), clamp(sp[1]+pad,DMIN,DMAX)]; }
  function pickTx(t){
    if (selectedOrf!=null) pickOrf(null);   // acting on a transcript ALWAYS escapes an ORF selection (single-writer)
    state[t].open=!state[t].open; if(!state[t].open) state[t].sharedOpen=false;
    selectedTx = state[t].open ? t : null;
    rebuild();   // expanding changes only the stack HEIGHT (not W) -> view is stable; no re-fit needed
    // ONE ACTION, ONE EFFECT: expanding a row is DISCLOSURE, not NAVIGATION. Do NOT change the
    // zoom/pan/locus — the 8 isoforms share ONE genomic axis so they can be compared, and the ORFs
    // render at their real coordinates within the current view. Navigation is deliberate: the locus
    // box / +- / scroll-pinch / ruler-brush / Reset. Name updates; coord stays = view.
    setNow(state[t].open ? ("MYCN-"+t) : "MYCN locus", null);
    renderFeatureList();
  }

  // ================= SINGLE-WRITER selection (FIX 1) =================
  // pickOrf(id|null) is the ONLY place `selectedOrf` is assigned in the whole file, so the variable and
  // the DOM can never diverge (that divergence — a control closing the panel without clearing the
  // variable — was what made the same-ORF toggle drift out of phase).
  // The Details panel PUSHES on desktop/tablet (in-flow column) and OVERLAYS as a full-width sheet on
  // phone. Opening/closing it is a WIDTH change, so we rebuild at the new width and re-apply the view —
  // the same "capture view -> rebuild -> applyRange" machinery a window resize uses. The genomic window
  // is preserved (the plot just gets narrower/wider); the view is NOT panned or zoomed (FIX 2).
  function pickOrf(id){
    // FIX 1(a): normalise the id type at the single boundary. Measurement showed every call site already
    // passes a String ORF id, but coercing here GUARANTEES the strict `id===selectedOrf` toggle test always
    // compares like types (never string-vs-number) regardless of any future datum-typed caller. Keys in
    // META.orfs are strings, so lookups are unaffected. (Never loose ==, per the stop-gate.)
    id = (id==null) ? null : String(id);
    if (id!=null && id===selectedOrf) id=null;   // clicking the SAME (selected) ORF toggles it off
    selectedOrf = id;                            // <-- the ONLY assignment to selectedOrf anywhere
    var target = viewG.slice();                  // default: keep the current genomic window across the resize
    if (id==null){
      railHasContent=false;
      $("rail-body").innerHTML='<div class="rail-empty">Select an ORF for its protein detail, or click an exon for its coordinates.</div>';
      openRight(false);
      setNow(selectedTx?("MYCN-"+selectedTx):"MYCN locus", selectedTx?txSpan(selectedTx):null);
    } else {
      var o=META.orfs[id], sp=o.span;
      $("rail-body").innerHTML = renderCard(id,o); wireCard();
      railHasContent=true;
      // FIX 4: expand the carrier the view scrolls to, so the selected ORF is actually VISIBLE on the stack
      // (e.g. after Collapse all) instead of hidden inside a collapsed transcript. Completes the action.
      var fc = TX_ORDER.filter(function(tt){ return carriers(id).indexOf(tt)>=0; })[0];
      if (fc) state[fc].open = true;
      var visible = (sp[1] >= viewG[0] && sp[0] <= viewG[1]);
      if (!visible) target = paddedSpan(sp);     // ORF off-screen -> bring it into view
      setNow("ORF"+id+(NAME[id]?" · "+NAME[id]:""), visible ? null : sp);
      openRight(true);                           // desktop: viewer narrows; phone: full-screen sheet
    }
    rebuild();                                   // recompute W at the new panel width; apply highlight/dim
    applyRange(target[0], target[1], false);     // re-fit the view at the new width (window preserved, or ORF span)
    scrollCarrierIntoView(id);
    markFeatureActive();
  }
  function showExon(t,exon,ei){
    var cds=META.transcripts[t].cds||[], isCds=false, phase=null;
    cds.forEach(function(c,ci){ if(!(c[1]<exon[0]||c[0]>exon[1])){ isCds=true; if(phase===null) phase=CDS_PHASE[t][ci]; } });
    var len=exon[1]-exon[0]+1, riHit=(RETAINED[t]||[]).some(function(iv){ return !(iv[1]<exon[0]||iv[0]>exon[1]); });
    var h='<div class="exon-detail"><h2>Exon detail</h2>';
    h+='<p class="exon-title">MYCN-'+esc(t)+' · exon '+(ei+1)+'</p><table class="props">';
    h+='<tr><td class="k">Coordinates</td><td class="v">chr2:'+fmt(exon[0])+'–'+fmt(exon[1])+'</td></tr>';
    h+='<tr><td class="k">Length</td><td class="v">'+fmt(len)+' bp</td></tr>';
    h+='<tr><td class="k">Content</td><td class="v">'+(isCds?"coding (CDS)":"non-coding (UTR)")+'</td></tr>';
    if (isCds) h+='<tr><td class="k">CDS phase (+ strand)</td><td class="v">'+phase+'</td></tr>';
    if (riHit) h+='<tr><td class="k">Retained intron</td><td class="v">yes ('+esc(t)+')</td></tr>';
    h+='<tr><td class="k">Strand</td><td class="v">'+esc(META.transcripts[t].strand)+'</td></tr></table>';
    h+='<p class="rail-empty">Click an ORF (expand the transcript) to see its protein detail.</p></div>';
    $("rail-body").innerHTML=h; railHasContent=true; openRight(true);
    setNow("MYCN-"+t+" exon "+(ei+1), exon);
    // FIX 2: opening Details is a WIDTH change on desktop/tablet — re-render at the new width NOW (mirror
    // pickOrf's tail) instead of waiting for the next resize, preserving the genomic window (no pan/zoom).
    var target = viewG.slice(); rebuild(); applyRange(target[0], target[1], false);
  }

  // tooltips
  function tip(ev,html){ var el=$("tip"); el.innerHTML=html; el.hidden=false; moveTip(ev); }
  function moveTip(ev){ var el=$("tip"); if(el.hidden) return;   // T2.3: clamp to viewport so tips don't overflow edges
    var vw=window.innerWidth||1200, vh=window.innerHeight||800, tw=el.offsetWidth||180, th=el.offsetHeight||40;
    var lx=ev.clientX+14, ty=ev.clientY+14;
    if (lx+tw+8>vw) lx=Math.max(8, ev.clientX-tw-14);
    if (ty+th+8>vh) ty=Math.max(8, ev.clientY-th-14);
    el.style.left=lx+"px"; el.style.top=ty+"px"; }
  function hideTip(){ $("tip").hidden=true; }
  function exonTip(t,exon,ei,s){ var kind=s.type==="cds"?"CDS":(s.type==="ri"?"retained intron":"UTR");
    return "<b>MYCN-"+t+" exon "+(ei+1)+"</b> · "+kind+"<br><span class='dim'>chr2:"+fmt(exon[0])+"-"+fmt(exon[1])+" · "+fmt(exon[1]-exon[0]+1)+" bp</span>"; }
  function orfTip(id){ var o=META.orfs[id], nm=NAME[id]?(" · "+NAME[id]):"";
    return "<b>ORF"+id+nm+"</b> · "+esc(o.cls)+"<br><span class='dim'>"+o.aa_len+" aa · "+carrierCount(id)+" carrier"+(carrierCount(id)===1?"":"s")+" · click for detail</span>"; }

  // ================= LEFT drawer: browse list (dual lens + search) =================
  function setLens(l){
    if (l===lens) return;
    lens=l;
    $("lens-tx").classList.toggle("active",l==="tx"); $("lens-orf").classList.toggle("active",l==="orf");
    $("cue").firstChild && ($("cue").childNodes[0].nodeValue = l==="tx"
      ? "Click a transcript to see the proteins it encodes." : "Click an ORF to view it on the stack.");
    renderFeatureList();
    var sq=$("feature-search").value; if (sq) filterList(sq);
  }
  function renderFeatureList(){ if (lens==="orf") renderOrfList(); else renderTxList(); markFeatureActive(); }

  function renderTxList(){
    novelOpen=false;
    var host=$("feature-list"), h='<div class="fgroup">';
    TX_ORDER.forEach(function(t){
      var n=orfIdsOfTx(t).length, coding=txCoding(t);
      h+='<button class="trow" data-tx="'+t+'" title="MYCN-'+t+'">'+
         '<span class="tchev">'+(state[t].open?"▾":"▸")+'</span>'+
         '<span class="tname">MYCN-'+esc(t)+'</span>'+
         '<span class="tmeta">'+n+' ORFs · <span class="'+(coding==="non-coding"?"tnc":"")+'">'+coding+'</span></span></button>';
    });
    h+='</div>'; host.innerHTML=h;
    Array.prototype.forEach.call(host.querySelectorAll(".trow"),function(btn){
      btn.addEventListener("click",function(){ pickTx(btn.getAttribute("data-tx")); });
    });
  }
  function orfRowHTML(id){
    var o=META.orfs[id], nm=NAME[id]?(" · "+NAME[id]):"", dot=String(o.conserved).toUpperCase()==="YES"?"on":"off";
    return '<button class="frow" data-id="'+id+'" title="fly to ORF'+id+'">'+
      '<span class="sw" style="background:'+o.color+'"></span>'+
      '<span class="ftxt">ORF'+id+nm+' · '+esc(o.cls)+' · '+o.aa_len+' aa</span>'+
      '<span class="ntag">'+carrierCount(id)+'c</span>'+
      '<span class="cdot '+dot+'" title="'+(dot==="on"?"conservation-supported":"not conservation-supported (a weak discriminator for short uORFs)")+'"></span></button>';
  }
  function renderOrfList(){
    novelOpen=false;
    var host=$("feature-list"), h="";
    GROUPS.forEach(function(gr){
      var ids=orfsOfClass(gr.key); if(!ids.length) return;
      h+='<div class="fgroup"><div class="fghead">'+esc(gr.title)+'</div>'+ids.map(orfRowHTML).join("")+'</div>';
    });
    var novelIds=NOVEL.reduce(function(a,n){ return a.concat(orfsOfClass(n.key)); },[]);
    h+='<div class="fgroup novelwrap"><button id="novel-toggle" class="novel-toggle" aria-expanded="false">'+
       '<span class="chev">▸</span> Show '+novelIds.length+' novel ORFs</button><div id="novel-body" class="novel-body" hidden>';
    NOVEL.forEach(function(n){ var ids=orfsOfClass(n.key);
      h+='<div class="fsub"><div class="fghead sub">'+esc(n.title)+' ('+ids.length+')</div>'+ids.map(orfRowHTML).join("")+'</div>'; });
    h+='</div></div>'; host.innerHTML=h;
    Array.prototype.forEach.call(host.querySelectorAll(".frow"),function(btn){
      btn.addEventListener("click",function(){ pickOrf(btn.getAttribute("data-id")); });
    });
    var ntog=$("novel-toggle"); if(ntog) ntog.addEventListener("click",function(){ setNovel(!novelOpen); });   // T2.3 null-guard
  }
  function setNovel(open){ novelOpen=open; var b=$("novel-body"); if(!b) return; b.hidden=!open;
    var t=$("novel-toggle"); t.setAttribute("aria-expanded",open?"true":"false"); t.querySelector(".chev").textContent=open?"▾":"▸"; }
  function filterList(q){
    q=(q||"").trim().toLowerCase();
    var host=$("feature-list"), anyNovel=false;
    Array.prototype.forEach.call(host.querySelectorAll(".frow"),function(btn){
      var id=btn.getAttribute("data-id"), o=META.orfs[id];
      var hay=["ORF"+id,NAME[id]||"",o.cls,(o.carriers||[]).join(" "),
               // T2.2: include the DISPLAYED conservation wording so search matches the UI
               (String(o.conserved).toUpperCase()==="YES"?"conserved conservation-supported":"not conserved not conservation-supported")+" conservation",
               o.acc||""].join(" ").toLowerCase();
      var hit=!q||hay.indexOf(q)>=0; btn.style.display=hit?"":"none";
      if (hit&&q&&btn.closest(".novel-body")) anyNovel=true;
    });
    Array.prototype.forEach.call(host.querySelectorAll(".trow"),function(btn){
      var t=btn.getAttribute("data-tx"), hay=("MYCN-"+t+" "+txCoding(t)).toLowerCase();
      btn.style.display=(!q||hay.indexOf(q)>=0)?"":"none";
    });
    if (q&&anyNovel&&!novelOpen) setNovel(true);
  }
  function markFeatureActive(){
    Array.prototype.forEach.call(document.querySelectorAll(".frow.active,.trow.active"),function(b){ b.classList.remove("active"); });
    if (lens==="orf"&&selectedOrf){ var b=document.querySelector('.frow[data-id="'+selectedOrf+'"]'); if(b) b.classList.add("active"); }
    if (lens==="tx"&&selectedTx){ var t=document.querySelector('.trow[data-tx="'+selectedTx+'"]'); if(t) t.classList.add("active"); }
  }

  // ================= detail card (identical format to the igv viewer) =================
  function wireCard(){
    var body=$("rail-body");
    Array.prototype.forEach.call(body.querySelectorAll("button.copy"),function(btn){
      btn.addEventListener("click",function(){ copySeq(btn); });
    });
  }
  function copySeq(btn){
    var txt=btn.getAttribute("data-seq")||"", old=btn.getAttribute("data-label")||btn.textContent;
    btn.setAttribute("data-label",old);
    function ok(){ btn.textContent="copied"; btn.classList.add("ok"); setTimeout(function(){ btn.textContent=old; btn.classList.remove("ok"); },1100); }
    function fail(){ btn.textContent="copy failed"; btn.classList.add("err"); setTimeout(function(){ btn.textContent=old; btn.classList.remove("err"); },1400); }
    if (navigator.clipboard && navigator.clipboard.writeText){
      navigator.clipboard.writeText(txt).then(ok, function(){ if(!legacyCopy(txt)) fail(); else ok(); });
    } else { if (legacyCopy(txt)) ok(); else fail(); }
  }
  function legacyCopy(txt){
    try{ var ta=document.createElement("textarea"); ta.value=txt; ta.style.position="fixed"; ta.style.opacity="0";
      document.body.appendChild(ta); ta.select(); var okc=document.execCommand&&document.execCommand("copy");
      document.body.removeChild(ta); return !!okc; }catch(e){ return false; }
  }
  function badgesFor(id,o){
    var b=[], cons=String(o.conserved).toUpperCase()==="YES";
    // conservation as EVIDENCE, not a verdict: absence is a neutral informational chip, not a fail flag
    b.push('<span class="badge '+(cons?"cons":"noncons")+'">'+(cons?"conservation-supported":"not conservation-supported")+"</span>");
    if (id==="9")  b.push('<span class="badge gold">known functional uORF · MYCNOT</span>');
    if (id==="10") b.push('<span class="badge gold">known functional uORF · MUSEP</span>');
    // T2.1: aa counts derived from o.aa_len (not hardcoded 593/575)
    if (id==="3")  b.push('<span class="badge canon">canonical '+o.aa_len+' aa · two-executor confirmed</span>');
    if (id==="4")  b.push('<span class="badge cand">second candidate ('+o.aa_len+' aa)</span>');
    if (o.cls==="N-Myc-frame-ext"&&id!=="3"&&id!=="4") b.push('<span class="badge canon">predicted N-term extension</span>');
    // T2.1: the 206 intron-retention set = ORFs unique to 206 (derived from carriers, not hardcoded)
    if (specGroups("206").uniq.indexOf(id)>=0) b.push('<span class="badge ri">206 intron-retention</span>');
    return b.join("");
  }
  function propRow(k,v){ return (v===null||v===undefined||v==="")?"":'<tr><td class="k">'+esc(k)+'</td><td class="v">'+esc(v)+"</td></tr>"; }
  function codonStrip(nt,aa){ nt=String(nt||""); aa=String(aa||""); var out=[];
    for (var i=0;i<nt.length;i+=3){ var codon=nt.substr(i,3), idx=i/3, res=idx<aa.length?aa[idx]:"*";
      var cls=i===0?"start":(idx>=aa.length?"stop":""); out.push('<span class="codon '+cls+'"><b>'+esc(codon)+'</b><i>'+esc(res)+"</i></span>"); }
    return out.join(""); }
  function renderCard(id,o){
    var swatch='<span class="orf-swatch" style="background:'+esc(o.color)+'"></span>';
    var span=o.span?("chr2:"+o.span[0].toLocaleString()+"–"+o.span[1].toLocaleString()):"", carr=(o.carriers||[]).join(", ");
    var cons=String(o.conserved).toUpperCase()==="YES";
    var h='<h2>ORF detail</h2>';
    h+='<p class="orf-id">'+swatch+"ORF "+esc(id)+(NAME[id]?' <span class="orf-alias">'+esc(NAME[id])+'</span>':"")+"</p>";
    h+='<p class="orf-class">'+esc(o.cls)+" · "+esc(o.location||"")+"</p>";
    // FIX 2: scope the uniqueness claim (unique = only within the 8 transcripts shown)
    if (carrierCount(id)===1) h+='<p class="scope-note">Unique among the 8 transcripts shown.</p>';
    h+='<div class="badges">'+badgesFor(id,o)+"</div>";
    // FIX 1: conservation caveat next to the badge, shown for conservation-negative ORFs (the honest case)
    if (!cons) h+='<p class="cons-note">Conservation is a weak discriminator for short uORFs: the known functional uORFs MYCNOT and MUSEP are both conservation-negative.</p>';
    h+="<table class=\"props\">";
    h+=propRow("Length",(o.aa_len!=null?o.aa_len+" aa":""));
    h+=propRow("Genomic span",span); h+=propRow("Carrier transcripts",carr);
    h+=propRow("Kozak (−3/+4)",o.kozak); h+=propRow("PhyloCSF",o.phylocsf); h+=propRow("phyloP (100-way)",o.phylop);
    h+=propRow("In Ribo-seq catalogue",o.catalogue); h+=propRow("MW (Da)",o.mw); h+=propRow("pI",o.pI);
    h+=propRow("GRAVY",o.gravy); h+=propRow("Net charge (pH 7)",o.netq); h+=propRow("Instability index",o.instab);
    h+="</table>";
    h+='<div class="seqblock"><div class="seqhead"><span class="lbl">codons</span></div><div class="codon-wrap">'+codonStrip(o.nt,o.aa)+"</div>";
    h+='<div class="codon-legend"><span class="s">start (ATG→M)</span><span class="e">stop</span></div></div>';
    h+='<div class="seqblock"><div class="seqhead"><span class="lbl">amino-acid sequence ('+(o.aa?o.aa.length:0)+' aa)</span>'+
       '<button class="copy" data-seq="'+esc(o.aa)+'">copy</button></div><div class="seq">'+esc(o.aa)+"</div></div>";
    h+='<div class="seqblock"><div class="seqhead"><span class="lbl">coding nt (ATG→stop, '+(o.nt?o.nt.length:0)+' nt)</span>'+
       '<button class="copy" data-seq="'+esc(o.nt)+'">copy</button></div><div class="seq">'+esc(o.nt)+"</div></div>";
    h+='<div class="provenance">Accession: <span class="acc">'+esc(o.acc||"—")+"</span><br>"+esc(o.source||"")+"</div>";
    return h;
  }

  // ================= export =================
  var SVGNS="http://www.w3.org/2000/svg";
  // T1.4: compose ideogram + ruler/reference + stack into ONE figure (stacked vertically, shared x),
  // preserving the current view/zoom/expansion. Returns {node,width,height}.
  function composeExport(){
    var GAP=6, totalH=IDEO_H+GAP+HDR_H+GAP+H;
    var out=document.createElementNS(SVGNS,"svg");   // namespace implicit via createElementNS (don't re-add xmlns -> would duplicate)
    out.setAttribute("width",W); out.setAttribute("height",totalH); out.setAttribute("viewBox","0 0 "+W+" "+totalH);
    var style=document.createElementNS(SVGNS,"style"); style.textContent="text{font-family:"+FONT+"}"; out.appendChild(style);
    var bg=document.createElementNS(SVGNS,"rect");
    bg.setAttribute("x",0); bg.setAttribute("y",0); bg.setAttribute("width",W); bg.setAttribute("height",totalH); bg.setAttribute("fill","#ffffff");
    out.appendChild(bg);
    // each source svg's children are copied into a <g> translated to its band (clipPaths move with the group)
    [[ "ideo",0 ],[ "hdr",IDEO_H+GAP ],[ "stack",IDEO_H+GAP+HDR_H+GAP ]].forEach(function(pair){
      var src=$(pair[0]); if(!src) return;
      var g=document.createElementNS(SVGNS,"g"); g.setAttribute("transform","translate(0,"+pair[1]+")");
      Array.prototype.forEach.call(src.childNodes,function(n){ g.appendChild(n.cloneNode(true)); });
      out.appendChild(g);
    });
    return { node:out, width:W, height:totalH };
  }
  function serializeExport(){
    var c=composeExport();
    return { xml:'<?xml version="1.0" encoding="UTF-8"?>\n'+new XMLSerializer().serializeToString(c.node), width:c.width, height:c.height };
  }
  function triggerDownload(blob,name){
    var url=URL.createObjectURL(blob), a=document.createElement("a");
    a.href=url; a.download=name; document.body.appendChild(a); a.click(); document.body.removeChild(a);
    setTimeout(function(){ URL.revokeObjectURL(url); },600);
  }
  // T1.3: UTF-8 safe base64 without deprecated unescape()
  function b64utf8(s){ return btoa(encodeURIComponent(s).replace(/%([0-9A-F]{2})/g,function(_,p){ return String.fromCharCode("0x"+p); })); }
  function downloadSVG(){
    try{ triggerDownload(new Blob([serializeExport().xml],{type:"image/svg+xml;charset=utf-8"}),"MYCN_isoform_stack.svg"); }
    catch(e){ console.error(e); toast("SVG export failed",true); }
  }
  function downloadPNG(){
    var ex, scale=2.5;
    try{ ex=serializeExport(); }catch(e){ console.error(e); toast("PNG export failed",true); return; }
    var img=new Image();
    img.onload=function(){
      try{
        var canvas=document.createElement("canvas");
        canvas.width=Math.round(ex.width*scale); canvas.height=Math.round(ex.height*scale);
        var ctx=canvas.getContext&&canvas.getContext("2d");
        if(!ctx){ toast("PNG export failed (no canvas)",true); return; }
        ctx.fillStyle="#ffffff"; ctx.fillRect(0,0,canvas.width,canvas.height);
        ctx.setTransform(scale,0,0,scale,0,0); ctx.drawImage(img,0,0);
        if(!canvas.toBlob){ toast("PNG export failed (no toBlob)",true); return; }
        canvas.toBlob(function(blob){ if(blob) triggerDownload(blob,"MYCN_isoform_stack.png"); else toast("PNG export failed",true); },"image/png");
      }catch(e){ console.error(e); toast("PNG export failed",true); }
    };
    img.onerror=function(){ toast("PNG export failed",true); };
    try{ img.src="data:image/svg+xml;base64,"+b64utf8(ex.xml); }
    catch(e){ console.error(e); toast("PNG export failed",true); }
  }

  // ================= headless self-test hooks =================
  try { if (typeof globalThis!=="undefined") globalThis.__mstack = {
    setMETA:function(m,s,c){ META=m; LOCSEQ=s||LOCSEQ; CYTO=c||CYTO; if(META.palette){ for(var k in META.palette) if(PAL[k]!==undefined) PAL[k]=META.palette[k]; } TX_ORDER.forEach(function(t){ state[t]={open:false,sharedOpen:false}; }); deriveIntrons(); deriveGPos(); buildLegendData(); },
    acenSpan:function(){ return acenSpan(); }, bandAt:function(p){ return bandAt(p); }, STAIN:STAIN,
    verifyBaseMap:verifyBaseMap, baseAt:baseAt,
    specGroups:specGroups, exonSegments:exonSegments, orfIdsOfTx:orfIdsOfTx, orfsOfClass:orfsOfClass,
    carrierCount:carrierCount, txCoding:txCoding, txSpan:txSpan, renderCard:renderCard, badgesFor:badgesFor,
    codonGenomic:function(id,k){ var g=GPOS[id]; return [g[3*k],g[3*k+1],g[3*k+2]]; },
    flyTo:function(a,b){ return flyTo(a,b,false); }, view:function(){ return viewG.slice(); }, ppb:pxPerBase,
    zoomToRange:function(a,b){ return zoomToRange(a,b,false); }, brushEnd:function(px0,px1){ return onRulerBrushEnd({selection:[px0,px1]}); },
    xOf:function(p){ return x(p); }, invX:function(px){ return x.invert(px); }, plotL:plotL, plotR:plotR,
    rightOpen:function(){ return rightOpen; }, openRight:openRight, noAnim:function(v){ NOANIM=!!v; },
    zoomFilter:zoomFilter, scrollStackBy:scrollStackBy,
    expandAll:function(full){ TX_ORDER.forEach(function(t){ state[t].open=true; if(full) state[t].sharedOpen=true; }); rebuild(); renderFeatureList(); },
    pickOrf:function(id){ return pickOrf(id); }, selectedOrf:function(){ return selectedOrf; },
    serializeExport:function(){ return serializeExport(); }, dmycnTx:dmycnTx,
    GPOS:function(){ return GPOS; }, RETAINED:function(){ return RETAINED; },
    CANON_INTRONS:function(){ return CANON_INTRONS; }, LEGEND:function(){ return LEGEND; }, PAL:function(){ return PAL; },
    TX_ORDER:TX_ORDER, BASECOL:BASECOL, GROUPS:GROUPS, NOVEL:NOVEL
  }; } catch(e){}
})();
