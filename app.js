(function() {
    "use strict";
    var PAL = {
        nmyc: "#4A6FA5",
        dmycn: "#D98A3D",
        mycnot: "#C0504D",
        musep: "#B0568C",
        ri_line: "#2F7A57",
        ri_bg: "#C3E0D4",
        pred_ext: "#8FA9C9",
        utr: "#C7CFCC",
        term: "#8A9791",
        gold: "#C9A227",
        evid: "#2B2BB5",
        ink: "#000000"
    };
    function msDetected(o) {
        return false;
    }
    function stopOf(o) {
        return o && o.blocks && o.blocks.length ? o.blocks[o.blocks.length - 1][1] : null;
    }
    function msStopGroups() {
        return {};
    }
    function msGroupOf() {
        return null;
    }
    var BACKBONE = "#b9c2c0";
    var FONT = "Inter, system-ui, -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif";
    var MONO = "ui-monospace, 'SF Mono', 'Roboto Mono', Menlo, Consolas, monospace";
    var BASECOL = {
        A: "#33A02C",
        C: "#1F78B4",
        G: "#FF7F00",
        T: "#E31A1C",
        U: "#E31A1C",
        N: "#B0B0B0"
    };
    var TX_ORDER = [ "201", "202", "203", "204", "205", "206", "207", "XM" ];
    var NAME = {};
    var GROUPS = [ {
        key: "canonical-CDS",
        title: "Database-annotated CDS"
    } ];
    var NOVEL = [ {
        key: "novel-uORF",
        title: "Upstream ORF predictions"
    }, {
        key: "novel-other",
        title: "Other ORF predictions"
    } ];
    var LEGEND = [];
    var GUTTER = 200, PADR = 18;
    var TXH = 64, ORFH = 23, GRPH = 20, SHH = 44, EXPAD = 8;
    var REGH = 52;
    function rowH(r) {
        return r.kind === "region" ? REGH : ORFH;
    }
    var THICK = 15, THIN = 6, RIH = 9, ORF_THICK = 11;
    var DMIN = 15938e3, DMAX = 15947200;
    var CHR2LEN = 242193529;
    var LETTER_PXPB = 7, CELL_PXPB = .5;
    var HDR_H = 56, IDEO_H = 56;
    var STAIN = {
        gneg: "#ffffff",
        gpos25: "#c9cdd1",
        gpos50: "#9aa1a7",
        gpos75: "#6f767c",
        gpos100: "#4a5157",
        gvar: "#c2c9cf",
        stalk: "#9aa7b0",
        acen: "#b06a80"
    };
    var LAB = window.MYCNLab;
    var META = null, LOCSEQ = null, CYTO = null;
    var x0 = null, x = null, zt = null, zoom = null, rulerBrush = null, dragPrevY = null;
    var svg = null, hdr = null, ideo = null, gMain = null;
    var W = 900, H = 400;
    var ROWS = [];
    var lens = "tx", selectedOrf = null, selectedTx = null, novelOpen = false, railHasContent = false, rightOpen = false;
    var selectedRegion = null;
    var viewG = [ DMIN, DMAX ];
    var state = {};
    var CANON_INTRONS = [], RETAINED = {}, CDS_PHASE = {}, GPOS = {};
    function $(id) {
        return document.getElementById(id);
    }
    function d3s(id) {
        return d3.select("#" + id);
    }
    function esc(s) {
        return String(s).replace(/[&<>"]/g, function(c) {
            return {
                "&": "&amp;",
                "<": "&lt;",
                ">": "&gt;",
                '"': "&quot;"
            }[c];
        });
    }
    function fmt(n) {
        return Math.round(n).toLocaleString("en-US");
    }
    function clamp(v, a, b) {
        return v < a ? a : v > b ? b : v;
    }
    function toast(msg, isErr) {
        var el = $("loadmsg");
        if (!el) return;
        if (!msg) {
            el.classList.add("hide");
            return;
        }
        el.textContent = msg;
        el.classList.remove("hide");
        el.classList.toggle("err", !!isErr);
    }
    window.addEventListener("DOMContentLoaded", function() {
        Promise.all([ fetch("data/current/reference.json", {cache:"no-cache"}).then(function(r) {
            if (!r.ok) throw new Error("meta.json HTTP " + r.status);
            return r.json();
        }), fetch("mycn_locus_seq.json").then(function(r) {
            if (!r.ok) throw new Error("locus_seq HTTP " + r.status);
            return r.json();
        }), fetch("chr2_cytoband.json").then(function(r) {
            if (!r.ok) throw new Error("cytoband HTTP " + r.status);
            return r.json();
        }) ]).then(function(res) {
            META = res[0];
            LOCSEQ = res[1];
            CYTO = res[2];
            TX_ORDER = META.transcript_order.slice();
            NAME = {};
            validateMetadata();
            if (!verifyBaseMap()) throw new Error("base-map ATG spot-check FAILED at chr2:15,942,065");
            init();
            buildDataSummary();
            toast(null);
        }).catch(function(e) {
            console.error(e);
            toast("Failed to load MYCN data: " + e.message, true);
        });
    });
    function validateMetadata() {
        if (!LAB || !META || !META.transcripts || !META.orfs) throw new Error("Missing transcript metadata or laboratory aliases");
        var ids = {};
        TX_ORDER.forEach(function(t) {
            var tx = META.transcripts[t];
            if (!tx || !tx.transcript_id || tx.id !== tx.transcript_id || ids[tx.id] || !tx.source || !tx.gencode_release || !Array.isArray(tx.exons) || !tx.exons.length || !Array.isArray(tx.cds)) throw new Error("Invalid transcript record: " + t);
            ids[tx.id] = true;
            var len = 0, last = 0;
            tx.exons.forEach(function(e) {
                if (e.length !== 2 || !Number.isInteger(e[0]) || !Number.isInteger(e[1]) || e[0] > e[1] || e[0] <= last) throw new Error("Invalid exon: " + t);
                last = e[1];
                len += e[1] - e[0] + 1;
            });
            if (len !== tx.mlen) throw new Error("Transcript length mismatch: " + t);
        });
        Object.keys(META.orfs).forEach(function(id) {
            var o = META.orfs[id];
            if (!o.carriers || !o.carriers.length || o.carriers.some(function(t) {
                return TX_ORDER.indexOf(t) < 0;
            }) || !o.blocks || !o.blocks.length || typeof o.aa !== "string" || o.aa.length !== o.aa_len) throw new Error("Invalid ORF: " + id);
        });
    }
    function buildDataSummary() {
        var body = $("transcript-summary");
        body.textContent = "";
        TX_ORDER.forEach(function(t) {
            var tx = META.transcripts[t], row = document.createElement("tr");
            [ LAB.label(t), txAcc(t), txSrc(t), txRel(t), String(tx.mlen), tx.exons.map(function(e) {
                return e.join("–");
            }).join("; "), tx.cds.length ? tx.cds.map(function(e) {
                return e.join("–");
            }).join("; ") : "No annotated CDS" ].forEach(function(v) {
                var cell = document.createElement("td");
                cell.textContent = v;
                row.appendChild(cell);
            });
            body.appendChild(row);
        });
        $("reference-summary").textContent = "CURRENT · GENCODE " + META.gencode_release + " · " + META.assembly + " · Ensembl " + META.ensembl_release + " · RefSeq " + META.refseq_accession;
        $("reference-version").textContent = "Scientific review: " + META.review_date + " · Application " + META.application_version + " · Reference build " + META.reference_build;
        $("reference-banner").textContent = "CURRENT · GENCODE " + META.gencode_release + " + RefSeq · " + META.assembly;
        $("reference-links").innerHTML = '<a href="https://www.gencodegenes.org/human/release_' + META.gencode_release + '.html">Official GENCODE release</a> · <a href="https://www.ncbi.nlm.nih.gov/nuccore/' + esc(META.refseq_accession) + '">Current RefSeq record</a> · <a href="https://www.uniprot.org/uniprotkb/P04198/entry">Reviewed UniProt P04198</a>';
        $("about-open").addEventListener("click", function() {
            $("about-data").showModal();
        });
        $("about-close").addEventListener("click", function() {
            $("about-data").close();
        });
        $("history-open").addEventListener("click", function() {
            $("release-history").showModal();
        });
        $("history-close").addEventListener("click", function() {
            $("release-history").close();
        });
        $("download-data").addEventListener("click", function() {
            var rows = TX_ORDER.map(function(t) {
                return Object.assign({
                    canonical_name: LAB.canonical(t),
                    lab_alias: LAB.alias(t)
                }, META.transcripts[t]);
            });
            triggerDownload(new Blob([ JSON.stringify({
                application_version: META.application_version,
                reference_build: META.reference_build,
                assembly: META.assembly,
                source_manifest: META.official_source_manifest,
                transcripts: rows,
                orfs: META.orfs
            }, null, 2) ], {
                type: "application/json"
            }), "MYCN_current_reference.json");
        });
        $("legend-toggle").addEventListener("keydown", function(e) {
            if (e.key === " " || e.key === "Enter") {
                e.preventDefault();
                e.currentTarget.click();
            }
        });
        document.addEventListener("keydown", function(e) {
            if (e.key === "Escape") {
                hideTip();
                var dialog = $("about-data").open ? $("about-data") : $("release-history").open ? $("release-history") : null;
                if (dialog) { e.preventDefault(); dialog.close(); $(dialog.id === "about-data" ? "about-open" : "history-open").focus(); return; }
                if (!dialog) {
                    pickOrf(null);
                    $("about-open").focus();
                }
            }
        });
    }
    function transcriptProperties(t) {
        var tx = META.transcripts[t], h = "";
        [ [ "Canonical name", tx.display_name ], [ "Laboratory alias", LAB.alias(t) || "Not supplied" ], [ "Biotype", tx.biotype ], [ "Transcript length", tx.mlen + " nt" ], [ "Exon count", tx.exons.length ], [ "CDS status", tx.cds.length ? "Annotated CDS" : "No annotated CDS" ], [ "Protein", tx.protein_accession ? tx.protein_accession + " · " + tx.protein_length + " aa" : "No supplied translation" ], [ "Annotation tags", (tx.tags || []).join(", ") || "None supplied" ], [ "Transcript support", tx.transcript_support_level ], [ "Evidence status", tx.evidence_class ], [ "Reference build", META.reference_build ], [ "Prior accession", tx.history ? tx.history.transcript_id : "Not in prior viewer" ], [ "Prior transcript length", tx.history ? tx.history.mlen + " nt" : "Not in prior viewer" ], [ "Prior exon coordinates", tx.history ? tx.history.exons.map(function(e){return e.join("–");}).join("; ") : "Not in prior viewer" ] ].forEach(function(r) {
            h += propRow(r[0], r[1]);
        });
        return h;
    }
    function proteinDiagram(o) {
        if (!o.domains.length) return "<p>No sequence-matched domain annotation assigned to this protein. This is not evidence that domains are absent.</p>";
        var h = '<svg viewBox="0 0 300 44" role="img" aria-label="Protein domains; exact positions listed below"><line x1="10" y1="15" x2="290" y2="15" stroke="#5f6a66" stroke-width="4"/>';
        o.domains.filter(function(d) {
            return d.type === "Domain";
        }).forEach(function(d) {
            var x = 10 + 280 * (d.start - 1) / o.aa_len, w = 280 * (d.end - d.start + 1) / o.aa_len;
            h += '<rect x="' + x + '" y="7" width="' + w + '" height="16" fill="#4A6FA5"/><text x="10" y="39" font-size="11">' + esc(d.description) + "</text>";
        });
        h += "</svg><ul>";
        o.domains.forEach(function(d) {
            h += "<li>" + esc(d.type + ": " + d.description) + " · " + d.start + "–" + d.end + " aa · " + esc(d.source) + " · " + esc(d.evidence_class) + "</li>";
        });
        return h + "</ul>";
    }
    function wireDiagramKeyboard() {
        document.querySelectorAll("#stack .gutter-hit").forEach(function(n) {
            var t = n.parentNode.getAttribute("data-tx");
            n.setAttribute("tabindex", "0");
            n.setAttribute("role", "button");
            n.setAttribute("aria-label", LAB.accessible(t) + ", " + txAcc(t) + ", " + txSrc(t) + ", reference " + txRel(t));
            n.setAttribute("aria-expanded", String(state[t].open));
            n.onkeydown = function(e) {
                if (e.key === "Enter" || e.key === " ") {
                    e.preventDefault();
                    pickTx(t);
                    var next = document.querySelector('#stack .tx-row[data-tx="' + t + '"] .gutter-hit');
                    if (next) next.focus();
                }
            };
        });
    }
    function verifyBaseMap() {
        if (!LOCSEQ || LOCSEQ.start !== DMIN) return false;
        var i = 15942065 - LOCSEQ.start;
        return LOCSEQ.seq.substr(i, 3).toUpperCase() === "ATG";
    }
    function baseAt(pos) {
        var i = pos - LOCSEQ.start;
        return i >= 0 && i < LOCSEQ.seq.length ? LOCSEQ.seq.charAt(i).toUpperCase() : "N";
    }
    function init() {
        if (META.palette) {
            for (var k in META.palette) if (PAL[k] !== undefined) PAL[k] = META.palette[k];
        }
        TX_ORDER.forEach(function(t) {
            state[t] = {
                open: false,
                sharedOpen: false
            };
        });
        deriveIntrons();
        deriveGPos();
        buildLegendData();
        buildHeader();
        buildLegend();
        wireChrome();
        svg = d3s("stack");
        hdr = d3s("hdr");
        ideo = d3s("ideo");
        zt = d3.zoomIdentity;
        setupZoom();
        setupTouch();
        renderFeatureList();
        wireLegendToggle();
        $("rail-body").innerHTML = renderLanding();
        wireCard();
        setLeft(!isPhone());
        openRight(false);
        rebuild();
        window.addEventListener("resize", onResizeDebounced);
    }
    function wireLegendToggle() {
        var lg = document.querySelector(".legend"), tt = document.querySelector(".legend-title");
        if (!lg || !tt) return;
        if (isPhone()) lg.classList.add("collapsed");
        tt.setAttribute("aria-expanded", String(!lg.classList.contains("collapsed")));
        tt.addEventListener("click", function() {
            lg.classList.toggle("collapsed");
            tt.setAttribute("aria-expanded", String(!lg.classList.contains("collapsed")));
        });
    }
    function carriers(id) {
        return META.orfs[id].carriers;
    }
    function carrierCount(id) {
        return META.orfs[id].carriers.length;
    }
    function orfIdsOfTx(t) {
        return Object.keys(META.orfs).filter(function(id) {
            return carriers(id).indexOf(t) >= 0;
        });
    }
    function orfsOfClass(cls) {
        return Object.keys(META.orfs).filter(function(id) {
            return META.orfs[id].cls === cls;
        }).sort(function(a, b) {
            return +a - +b;
        });
    }
    function txAcc(t) {
        var x = META.transcripts[t] || {};
        return x.transcript_id || "[no accession on record]";
    }
    function txSrc(t) {
        var x = META.transcripts[t] || {};
        return x.source || "[source not on record]";
    }
    function txRel(t) {
        var x = META.transcripts[t] || {};
        return x.reference_release || x.gencode_release || "[release not on record]";
    }
    function txOffV49() {
        return false;
    }
    function txMane(t) {
        return (META.transcripts[t].mane || []).indexOf("MANE_Select") >= 0;
    }
    function txTagShort(t) {
        return t === "XM" ? " · RefSeq" : " · " + txRel(t);
    }
    function txTag(t) {
        return " · " + txSrc(t);
    }
    function eventCountOfTx(t) {
        var ids = orfIdsOfTx(t), seen = {}, n = 0;
        ids.forEach(function(id) {
            var rid = regionOf(id);
            if (rid && regionCollapses(rid)) {
                if (seen[rid]) return;
                seen[rid] = 1;
            }
            n++;
        });
        return n;
    }
    function txCoding(t) {
        var c = META.transcripts[t].cds;
        return c && c.length ? "coding" : "non-coding";
    }
    function txSpan(t) {
        var e = META.transcripts[t].exons;
        return [ e[0][0], e[e.length - 1][1] ];
    }
    function deriveIntrons() {
        var ex = META.transcripts["201"].exons;
        CANON_INTRONS = [];
        for (var i = 0; i < ex.length - 1; i++) CANON_INTRONS.push([ ex[i][1] + 1, ex[i + 1][0] - 1 ]);
        TX_ORDER.forEach(function(t) {
            RETAINED[t] = [];
            META.transcripts[t].exons.forEach(function(e) {
                CANON_INTRONS.forEach(function(iv) {
                    if (e[0] < iv[0] && e[1] > iv[1]) RETAINED[t].push(iv.slice());
                });
            });
            var cds = META.transcripts[t].cds || [], cum = 0;
            CDS_PHASE[t] = cds.map(function(c) {
                var ph = (3 - cum % 3) % 3;
                cum += c[1] - c[0] + 1;
                return ph;
            });
        });
    }
    function deriveGPos() {
        Object.keys(META.orfs).forEach(function(id) {
            var g = [];
            META.orfs[id].blocks.forEach(function(b) {
                for (var p = b[0]; p <= b[1]; p++) g.push(p);
            });
            GPOS[id] = g;
        });
    }
    function specGroups(t) {
        var mem = orfIdsOfTx(t), uniq = [], shared = [];
        mem.forEach(function(id) {
            if (carrierCount(id) === 1) uniq.push(id); else shared.push(id);
        });
        var byC = function(a, b) {
            return carrierCount(a) - carrierCount(b) || +a - +b;
        };
        uniq.sort(byC);
        shared.sort(byC);
        return {
            uniq: uniq,
            shared: shared
        };
    }
    function exonSegments(t, exon) {
        var cds = META.transcripts[t].cds || [], ri = RETAINED[t] || [];
        var cuts = {};
        cuts[exon[0]] = 1;
        cuts[exon[1] + 1] = 1;
        function addCuts(ivs) {
            ivs.forEach(function(iv) {
                var a = Math.max(iv[0], exon[0]), b = Math.min(iv[1], exon[1]);
                if (a <= b) {
                    cuts[a] = 1;
                    cuts[b + 1] = 1;
                }
            });
        }
        addCuts(cds);
        addCuts(ri);
        var pts = Object.keys(cuts).map(Number).sort(function(a, b) {
            return a - b;
        }), segs = [];
        for (var i = 0; i < pts.length - 1; i++) {
            var a = pts[i], b = pts[i + 1] - 1;
            if (b < a) continue;
            var mid = (a + b) / 2, type = "utr";
            if (inAny(cds, mid)) type = "cds"; else if (inAny(ri, mid)) type = "ri";
            segs.push({
                a: a,
                b: b,
                type: type
            });
        }
        return segs;
    }
    function inAny(ivs, p) {
        for (var i = 0; i < ivs.length; i++) if (p >= ivs[i][0] && p <= ivs[i][1]) return true;
        return false;
    }
    function buildHeader() {
        var nT = TX_ORDER.length, nO = Object.keys(META.orfs).length;
        $("tb-head").textContent = "MYCN locus";
        $("tb-sub").textContent = "— " + nT + " transcripts · " + nO + " ORFs";
    }
    function buildLegendData() {
        LEGEND = [ [ PAL.nmyc, "Annotated coding sequence (CDS)", false ], [ PAL.dmycn, "Alternative annotated CDS", false ], [ PAL.utr, "Untranslated or no-CDS region", true ], [ PAL.ri_line, "Intron-spanning exon relative to MYCN-201", true ], [ "#5B7D71", "Upstream ORF prediction", false ], [ "#657D96", "Other ORF prediction", false ] ];
    }
    function buildLegend() {
        $("legend-list").innerHTML = LEGEND.map(function(p) {
            return '<li><i style="background:' + p[0] + '"></i>' + esc(p[1]) + "</li>";
        }).join("");
    }
    function wireChrome() {
        var cx = $("cue-dismiss");
        if (cx) cx.addEventListener("click", function() {
            $("cue").hidden = true;
        });
        $("expand-all").addEventListener("click", function() {
            TX_ORDER.forEach(function(t) {
                state[t].open = true;
                state[t].sharedOpen = true;
            });
            rebuild();
            renderFeatureList();
        });
        $("collapse-all").addEventListener("click", function() {
            TX_ORDER.forEach(function(t) {
                state[t].open = false;
                state[t].sharedOpen = false;
            });
            selectedTx = null;
            pickOrf(null);
            renderFeatureList();
        });
        $("download-svg").addEventListener("click", downloadSVG);
        $("download-png").addEventListener("click", downloadPNG);
        $("locus-form").addEventListener("submit", function(ev) {
            ev.preventDefault();
            gotoLocusInput($("locus-input").value);
        });
        $("zoom-in").addEventListener("click", function() {
            zoomButton(1.6);
        });
        $("zoom-out").addEventListener("click", function() {
            zoomButton(1 / 1.6);
        });
        $("reset-view").addEventListener("click", resetView);
        $("lens-tx").addEventListener("click", function() {
            setLens("tx");
        });
        $("lens-orf").addEventListener("click", function() {
            setLens("orf");
        });
        $("feature-search").addEventListener("input", function() {
            filterList(this.value);
        });
        $("features-strip").addEventListener("click", function() {
            var v = viewG.slice();
            setLeft(true);
            rebuild();
            applyRange(v[0], v[1], false);
        });
        $("features-collapse").addEventListener("click", function() {
            var v = viewG.slice();
            setLeft(false);
            rebuild();
            applyRange(v[0], v[1], false);
        });
        $("details-strip").addEventListener("click", function() {
            if (railHasContent) openRight(true);
        });
        $("details-close").addEventListener("click", function() {
            pickOrf(null);
        });
    }
    function setLeft(open) {
        if (open && isPhone() && selectedOrf != null) pickOrf(null);
        if (open && isPhone() && rightOpen) openRight(false);
        $("features").hidden = !open;
        $("features-strip").hidden = open;
        syncDrawerAccessibility();
    }
    function openRight(open) {
        rightOpen = !!open;
        if (rightOpen && isPhone()) setLeft(false);
        $("details").hidden = !rightOpen;
        $("details-strip").hidden = rightOpen;
        syncDrawerAccessibility();
    }
    function syncDrawerAccessibility() {
        var phone = isPhone(), left = !$("features").hidden, right = !$("details").hidden, overlay = phone && (left || right);
        [ "topbar", "viewer", "features-strip", "details-strip" ].forEach(function(id) {
            $(id).inert = overlay;
        });
        [ [ "features", left, "Transcript and ORF features" ], [ "details", right, "Feature details" ] ].forEach(function(v) {
            var el = $(v[0]);
            el.setAttribute("role", phone && v[1] ? "dialog" : "complementary");
            el.setAttribute("aria-label", v[2]);
            if (phone && v[1]) el.setAttribute("aria-modal", "true"); else el.removeAttribute("aria-modal");
        });
        var active = document.activeElement;
        if (overlay && !$(left ? "features" : "details").contains(active)) $(left ? "feature-search" : "details-close").focus(); else if (!overlay && active && active.closest && active.closest("aside[hidden]")) $(phone ? "features-strip" : "about-open").focus();
    }
    function plotL() {
        return GUTTER;
    }
    function plotR() {
        return W - PADR;
    }
    function setupZoom() {
        zoom = d3.zoom().clickDistance(6).touchable(function() {
            return false;
        }).filter(zoomFilter).on("start", function(ev) {
            dragPrevY = ev.sourceEvent && ev.sourceEvent.clientY != null ? ev.sourceEvent.clientY : null;
        }).on("zoom", onZoom).on("end", function() {
            dragPrevY = null;
            try {
                svg.property("__zoom", zt);
            } catch (e) {}
        });
        svg.call(zoom);
    }
    var TCH = null, TAP_SLOP = 8;
    function svgScaleX() {
        var r = svg.node().getBoundingClientRect();
        return r.width ? W / r.width : 1;
    }
    function svgUserX(clientX) {
        var r = svg.node().getBoundingClientRect();
        return (clientX - r.left) * (r.width ? W / r.width : 1);
    }
    function setupTouch() {
        var node = svg.node();
        node.addEventListener("touchstart", tStart, {
            passive: false
        });
        node.addEventListener("touchmove", tMove, {
            passive: false
        });
        node.addEventListener("touchend", tEnd, {
            passive: false
        });
        node.addEventListener("touchcancel", tEnd, {
            passive: false
        });
    }
    function tStart(ev) {
        if (ev.touches.length === 1) {
            var t = ev.touches[0], sc = $("stack-scroll");
            TCH = {
                mode: "one",
                x0: t.clientX,
                y0: t.clientY,
                view0: viewG.slice(),
                st0: sc ? sc.scrollTop : 0,
                moved: false
            };
        } else if (ev.touches.length >= 2) {
            var a = ev.touches[0], b = ev.touches[1], mid = (a.clientX + b.clientX) / 2;
            TCH = {
                mode: "pinch",
                dist0: Math.max(1, Math.hypot(a.clientX - b.clientX, a.clientY - b.clientY)),
                view0: viewG.slice(),
                anchorG: x.invert(svgUserX(mid)),
                moved: true
            };
            ev.preventDefault();
        }
    }
    function tMove(ev) {
        if (!TCH) return;
        if (TCH.mode === "one" && ev.touches.length === 1) {
            var t = ev.touches[0], dx = t.clientX - TCH.x0, dy = t.clientY - TCH.y0;
            if (!TCH.moved && Math.hypot(dx, dy) < TAP_SLOP) return;
            TCH.moved = true;
            ev.preventDefault();
            var wbp = TCH.view0[1] - TCH.view0[0];
            var bpPerUser = wbp / Math.max(1, plotR() - plotL());
            var a = TCH.view0[0] - dx * svgScaleX() * bpPerUser;
            applyRange(a, a + wbp, false);
            var sc = $("stack-scroll");
            if (sc) sc.scrollTop = Math.max(0, TCH.st0 - dy);
        } else if (TCH.mode === "pinch" && ev.touches.length >= 2) {
            ev.preventDefault();
            var a2 = ev.touches[0], b2 = ev.touches[1], mid = (a2.clientX + b2.clientX) / 2;
            var dist = Math.max(1, Math.hypot(a2.clientX - b2.clientX, a2.clientY - b2.clientY));
            var wbp = clamp((TCH.view0[1] - TCH.view0[0]) / (dist / TCH.dist0), 20, DMAX - DMIN);
            var f = clamp((svgUserX(mid) - plotL()) / Math.max(1, plotR() - plotL()), 0, 1);
            var a = TCH.anchorG - f * wbp;
            applyRange(a, a + wbp, false);
        }
    }
    function tEnd(ev) {
        if (TCH && TCH.mode === "pinch" && ev.touches.length === 1) {
            var t = ev.touches[0], sc = $("stack-scroll");
            TCH = {
                mode: "one",
                x0: t.clientX,
                y0: t.clientY,
                view0: viewG.slice(),
                st0: sc ? sc.scrollTop : 0,
                moved: true
            };
            return;
        }
        if (ev.touches.length === 0) TCH = null;
    }
    function zoomFilter(ev) {
        if (ev.type === "wheel") return true;
        if (ev.type === "dblclick") return false;
        return !ev.button;
    }
    function scrollStackBy(dy) {
        var sc = $("stack-scroll");
        if (!sc) return 0;
        sc.scrollTop = Math.max(0, sc.scrollTop - dy);
        return sc.scrollTop;
    }
    function refreshZoomExtent() {
        var maxk = (DMAX - DMIN) / 40;
        zoom.extent([ [ plotL(), 0 ], [ plotR(), Math.max(1, H) ] ]).scaleExtent([ 1, Math.max(2, maxk) ]).translateExtent([ [ plotL(), -Infinity ], [ plotR(), Infinity ] ]);
        svg.call(zoom.transform, zt);
    }
    function onZoom(ev) {
        zt = d3.zoomIdentity.translate(ev.transform.x, 0).scale(ev.transform.k);
        var se = ev.sourceEvent;
        if (se && dragPrevY != null && se.clientY != null && (se.type === "mousemove" || se.type === "pointermove")) {
            var dy = se.clientY - dragPrevY;
            dragPrevY = se.clientY;
            if (dy) scrollStackBy(dy);
        }
        scheduleDraw();
    }
    var drawPending = false, RAF = typeof requestAnimationFrame === "function" ? requestAnimationFrame : function(f) {
        return setTimeout(f, 16);
    };
    function scheduleDraw() {
        if (NOANIM) {
            draw();
            return;
        }
        if (drawPending) return;
        drawPending = true;
        RAF(function() {
            drawPending = false;
            draw();
        });
    }
    function zoomButton(factor) {
        var c = (plotL() + plotR()) / 2;
        svg.interrupt();
        if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) svg.call(zoom.scaleBy, factor, [ c, 0 ]); else svg.transition().duration(180).call(zoom.scaleBy, factor, [ c, 0 ]);
    }
    function flyTo(g0, g1, animate) {
        var pad = Math.max(30, (g1 - g0) * .12);
        applyRange(clamp(g0 - pad, DMIN, DMAX), clamp(g1 + pad, DMIN, DMAX), animate);
    }
    function zoomToRange(g0, g1, animate) {
        if (g1 < g0) {
            var s = g0;
            g0 = g1;
            g1 = s;
        }
        applyRange(clamp(g0, DMIN, DMAX), clamp(g1, DMIN, DMAX), animate);
    }
    function applyRange(a, b, animate) {
        if (b - a < 20) {
            var m = (a + b) / 2;
            a = clamp(m - 10, DMIN, DMAX);
            b = clamp(m + 10, DMIN, DMAX);
        }
        var k = clamp((plotR() - plotL()) / (x0(b) - x0(a)), 1, zoom.scaleExtent()[1]);
        var tx = plotL() - k * x0(a);
        var t = d3.zoomIdentity.translate(tx, 0).scale(k);
        svg.interrupt();
        if (animate && !NOANIM && !window.matchMedia("(prefers-reduced-motion: reduce)").matches) svg.transition().duration(360).call(zoom.transform, t); else svg.call(zoom.transform, t);
    }
    var NOANIM = false;
    function resetView() {
        selectedTx = null;
        pickOrf(null);
        flyTo(DMIN, DMAX, true);
        setNow("MYCN locus", null);
        renderFeatureList();
    }
    function gotoLocusInput(str) {
        var m = String(str).replace(/,/g, "").match(/(?:chr2[:\s])?\s*(\d+)\s*[-–:]\s*(\d+)/i);
        if (!m) {
            var el = $("locus-input");
            if (el) {
                el.classList.add("bad");
                setTimeout(function() {
                    el.classList.remove("bad");
                }, 900);
            }
            return;
        }
        var a = +m[1], b = +m[2];
        if (a > b) {
            var t = a;
            a = b;
            b = t;
        }
        zoomToRange(clamp(a, DMIN, DMAX), clamp(b, DMIN, DMAX), true);
    }
    function vpW() {
        return typeof window !== "undefined" && window.innerWidth ? window.innerWidth : 1024;
    }
    function isPhone() {
        return vpW() < 640;
    }
    function isTablet() {
        var w = vpW();
        return w >= 640 && w < 1024;
    }
    function computeGutter() {
        if (isPhone()) {
            GUTTER = 132;
            PADR = 8;
        } else if (isTablet()) {
            GUTTER = 150;
            PADR = 12;
        } else {
            GUTTER = 200;
            PADR = 18;
        }
    }
    function computeW() {
        computeGutter();
        var box = $("stack-scroll"), cw = box ? box.clientWidth : 900;
        var floor = isPhone() ? 300 : isTablet() ? 480 : 720;
        W = Math.max(floor, cw - 2);
    }
    var SHORT_LBL_MAX = 130;
    function shortLabels() {
        return GUTTER < SHORT_LBL_MAX;
    }
    function layout() {
        try {
            TXH = txRowH();
            ORFH = orfRowH();
        } catch (e) {}
        var rows = [], y = 8;
        TX_ORDER.forEach(function(t) {
            rows.push({
                kind: "tx",
                t: t,
                y: y
            });
            y += TXH;
            if (state[t].open) {
                var g = specGroups(t);
                y = pushGroup(rows, "Unique ORFs", g.uniq, true, t, y);
                if (g.shared.length) {
                    rows.push({
                        kind: "shared",
                        t: t,
                        y: y,
                        n: g.shared.length
                    });
                    y += SHH;
                    if (state[t].sharedOpen) eventsOf(g.shared).forEach(function(e) {
                        if (e.kind === "region") {
                            rows.push({
                                kind: "region",
                                t: t,
                                rid: e.rid,
                                ids: e.ids,
                                uniq: false,
                                y: y
                            });
                            y += REGH;
                        } else {
                            rows.push({
                                kind: "orf",
                                t: t,
                                id: e.id,
                                uniq: false,
                                y: y
                            });
                            y += ORFH;
                        }
                    });
                }
                y += EXPAD;
            }
        });
        H = y + 12;
        return rows;
    }
    function pushGroup(rows, label, ids, uniq, t, y) {
        if (!ids.length) return y;
        rows.push({
            kind: "grp",
            label: label,
            uniq: uniq,
            y: y
        });
        y += GRPH;
        eventsOf(ids).forEach(function(e) {
            if (e.kind === "region") {
                rows.push({
                    kind: "region",
                    t: t,
                    rid: e.rid,
                    ids: e.ids,
                    uniq: uniq,
                    y: y
                });
                y += REGH;
            } else {
                rows.push({
                    kind: "orf",
                    t: t,
                    id: e.id,
                    uniq: uniq,
                    y: y
                });
                y += ORFH;
            }
        });
        return y;
    }
    function rebuild() {
        computeW();
        x0 = d3.scaleLinear().domain([ DMIN, DMAX ]).range([ plotL(), plotR() ]);
        ROWS = layout();
        svg.attr("width", W).attr("height", H).attr("viewBox", "0 0 " + W + " " + H).attr("font-family", FONT);
        hdr.attr("width", W).attr("height", HDR_H).attr("viewBox", "0 0 " + W + " " + HDR_H).attr("font-family", FONT);
        ideo.attr("width", W).attr("height", IDEO_H).attr("viewBox", "0 0 " + W + " " + IDEO_H).attr("font-family", FONT);
        refreshZoomExtent();
        drawIdeogramBar();
        draw();
    }
    function onResize() {
        syncDrawerAccessibility();
        var v = viewG.slice();
        rebuild();
        applyRange(v[0], v[1], false);
    }
    var resizeTimer = null;
    function onResizeDebounced() {
        if (resizeTimer) clearTimeout(resizeTimer);
        resizeTimer = setTimeout(onResize, 120);
    }
    function draw() {
        x = zt.rescaleX(x0);
        viewG = [ clamp(x.invert(plotL()), DMIN, DMAX), clamp(x.invert(plotR()), DMIN, DMAX) ];
        renderHeader();
        renderStack();
        updateIdeoMarker();
        updateNow();
        wireDiagramKeyboard();
        $("locus-input").value = "chr2:" + fmt(viewG[0]) + "-" + fmt(viewG[1]);
    }
    function chevrons(g, xa, xb, mid, color) {
        for (var cx = xa + 26; cx < xb - 5; cx += 52) {
            g.append("path").attr("d", "M" + (cx - 3) + " " + (mid - 4) + " L" + (cx + 3) + " " + mid + " L" + (cx - 3) + " " + (mid + 4)).attr("fill", "none").attr("stroke", color).attr("stroke-width", 1.1);
        }
    }
    function pxPerBase() {
        return (plotR() - plotL()) / (viewG[1] - viewG[0]);
    }
    function ideoGeom() {
        var chrLen = CYTO && CYTO.length || CHR2LEN;
        var L = GUTTER, R = W - PADR, midY = 26, barH = 14;
        var ix = function(p) {
            return L + (R - L) * (clamp(p, 0, chrLen) / chrLen);
        };
        return {
            chrLen: chrLen,
            L: L,
            R: R,
            midY: midY,
            barH: barH,
            ix: ix
        };
    }
    function acenSpan() {
        var a = (CYTO.bands || []).filter(function(b) {
            return b.stain === "acen";
        });
        if (!a.length) return null;
        return {
            start: a[0].start,
            mid: a[0].end,
            end: a[a.length - 1].end
        };
    }
    function bandAt(pos) {
        var bs = CYTO && CYTO.bands || [];
        for (var i = 0; i < bs.length; i++) if (pos >= bs[i].start && pos < bs[i].end) return bs[i];
        return null;
    }
    function drawIdeogramBar() {
        ideo.selectAll("*").remove();
        if (!CYTO || !CYTO.bands) {
            return;
        }
        var G = ideoGeom(), ix = G.ix, midY = G.midY, barH = G.barH, top = midY - barH / 2, bot = midY + barH / 2;
        var acen = acenSpan();
        ideo.append("text").attr("x", 12).attr("y", midY + 4).attr("font-size", 11).attr("font-weight", 600).attr("fill", "#5f6a66").text("chr2");
        var defs = ideo.append("defs");
        defs.append("clipPath").attr("id", "ideo-clip").append("rect").attr("x", ix(0)).attr("y", top).attr("width", ix(G.chrLen) - ix(0)).attr("height", barH).attr("rx", barH / 2).attr("ry", barH / 2);
        var gb = ideo.append("g").attr("clip-path", "url(#ideo-clip)");
        CYTO.bands.forEach(function(b) {
            if (b.stain === "acen") return;
            var xa = ix(b.start), xb = ix(b.end);
            gb.append("rect").attr("x", xa).attr("y", top).attr("width", Math.max(.4, xb - xa)).attr("height", barH).attr("fill", STAIN[b.stain] || "#d9dee2").append("title").text("chr2 " + b.name + " (" + b.stain + ")");
        });
        if (acen) {
            var xs = ix(acen.start), xm = ix(acen.mid), xe = ix(acen.end);
            gb.append("path").attr("d", "M" + xs + " " + top + " L" + xm + " " + midY + " L" + xe + " " + top + " Z").attr("fill", STAIN.acen);
            gb.append("path").attr("d", "M" + xs + " " + bot + " L" + xm + " " + midY + " L" + xe + " " + bot + " Z").attr("fill", STAIN.acen);
        }
        ideo.append("rect").attr("x", ix(0)).attr("y", top).attr("width", ix(G.chrLen) - ix(0)).attr("height", barH).attr("rx", barH / 2).attr("ry", barH / 2).attr("fill", "none").attr("stroke", "#9aa4ab").attr("stroke-width", 1);
        ideo.append("path").attr("id", "ideo-mark-tip").attr("fill", "#c0392b");
        ideo.append("rect").attr("id", "ideo-mark").attr("y", top - 3).attr("height", barH + 6).attr("rx", 2).attr("fill", "none").attr("stroke", "#c0392b").attr("stroke-width", 1.7);
        ideo.append("text").attr("id", "ideo-band").attr("y", bot + 13).attr("font-size", 9.5).attr("font-weight", 600).attr("fill", "#a5322c").attr("text-anchor", "middle").text("");
        updateIdeoMarker();
    }
    function updateIdeoMarker() {
        if (!CYTO || !CYTO.bands) return;
        var G = ideoGeom(), ix = G.ix, midY = G.midY, barH = G.barH, top = midY - barH / 2;
        var mid = (viewG[0] + viewG[1]) / 2;
        var xa = ix(viewG[0]), xb = ix(viewG[1]), cx = ix(mid);
        var w = Math.max(6, xb - xa);
        var m = ideo.select("#ideo-mark");
        if (m.empty()) return;
        m.attr("x", cx - w / 2).attr("width", w);
        ideo.select("#ideo-mark-tip").attr("d", "M" + (cx - 4) + " " + (top - 6) + " L" + (cx + 4) + " " + (top - 6) + " L" + cx + " " + (top - 1) + " Z");
        var band = bandAt(mid), name = band ? band.name : "";
        var halfW = name.length * .62 * 9.5 / 2 + 3;
        var lx = clamp(cx, halfW + 2, W - halfW - 2);
        ideo.select("#ideo-band").attr("x", lx).text(name);
    }
    function renderHeader() {
        hdr.selectAll("*").remove();
        hdr.append("defs").append("clipPath").attr("id", "hdr-clip").append("rect").attr("x", plotL()).attr("y", 0).attr("width", Math.max(1, plotR() - plotL())).attr("height", HDR_H);
        var rulerY = 14, stripY = 24, stripH = 13, codonY = stripY + stripH + 2;
        hdr.append("rect").attr("x", 0).attr("y", 0).attr("width", plotL() - 1).attr("height", HDR_H).attr("fill", "#ffffff");
        var g = hdr.append("g").attr("clip-path", "url(#hdr-clip)");
        g.append("line").attr("x1", plotL()).attr("x2", plotR()).attr("y1", rulerY).attr("y2", rulerY).attr("stroke", "#d7dde1");
        var plotW = plotR() - plotL();
        var target = Math.max(2, Math.min(8, Math.round(plotW / 90)));
        var lastRight = -1e9;
        x.ticks(target).forEach(function(v) {
            var xp = x(v);
            if (xp < plotL() - 1 || xp > plotR() + 1) return;
            g.append("line").attr("x1", xp).attr("x2", xp).attr("y1", rulerY - 4).attr("y2", rulerY).attr("stroke", "#c4ccd1");
            var txt = g.append("text").attr("class", "ruler-lbl").attr("x", xp).attr("y", rulerY - 6).attr("text-anchor", "middle").attr("font-size", 9.5).attr("font-family", MONO).attr("fill", "#5f6a66").text((v / 1e3).toFixed(2) + " kb");
            var w = txt.node().getComputedTextLength ? txt.node().getComputedTextLength() : 0;
            if (w && xp - w / 2 < lastRight + 8) txt.remove(); else lastRight = xp + (w || 0) / 2;
        });
        hdr.append("text").attr("x", 12).attr("y", stripY + stripH - 3).attr("font-size", 9.5).attr("font-weight", 600).attr("fill", "#5f6a66").text("ref");
        renderReferenceStrip(hdr.append("g").attr("clip-path", "url(#hdr-clip)"), stripY, stripH);
        renderCodonTrack(hdr.append("g").attr("clip-path", "url(#hdr-clip)"), codonY);
        installRulerBrush(codonY + 13);
    }
    function installRulerBrush(yBottom) {
        rulerBrush = d3.brushX().extent([ [ plotL(), 1 ], [ plotR(), Math.max(20, yBottom) ] ]).on("end", onRulerBrushEnd);
        var gB = hdr.append("g").attr("class", "ruler-brush").call(rulerBrush);
        gB.select(".overlay").attr("cursor", "ew-resize");
    }
    function onRulerBrushEnd(ev) {
        if (!ev.selection) return;
        var px0 = ev.selection[0], px1 = ev.selection[1];
        hdr.select(".ruler-brush").call(rulerBrush.move, null);
        if (px1 - px0 < 4) return;
        zoomToRange(x.invert(px0), x.invert(px1), true);
    }
    function renderReferenceStrip(g, y, h) {
        var ppb = pxPerBase();
        var a = Math.max(DMIN, Math.floor(viewG[0])), b = Math.min(DMAX, Math.ceil(viewG[1]));
        if (ppb >= LETTER_PXPB) {
            for (var p = a; p <= b; p++) {
                var base = baseAt(p), xa = x(p), xb = x(p + 1), cx = (xa + xb) / 2;
                g.append("rect").attr("x", xa).attr("y", y).attr("width", Math.max(1, xb - xa)).attr("height", h).attr("fill", BASECOL[base] || BASECOL.N).attr("opacity", .2);
                g.append("text").attr("class", "base-let").attr("x", cx).attr("y", y + h - 2.5).attr("text-anchor", "middle").attr("font-size", Math.min(13, (xb - xa) * .9)).attr("font-family", MONO).attr("font-weight", 600).attr("fill", BASECOL[base] || BASECOL.N).text(base);
            }
        } else {
            var step = ppb >= CELL_PXPB ? 1 : Math.max(1, Math.ceil((b - a) / 800));
            for (var q = a; q <= b; q += step) {
                var bb = baseAt(q), x1 = x(q), x2 = x(q + step);
                g.append("rect").attr("x", x1).attr("y", y).attr("width", Math.max(.6, x2 - x1)).attr("height", h).attr("fill", BASECOL[bb] || BASECOL.N).attr("opacity", ppb >= CELL_PXPB ? .9 : .75);
            }
        }
    }
    function renderCodonTrack(g, y) {
        if (pxPerBase() < LETTER_PXPB || !selectedOrf) return;
        var id = selectedOrf, o = META.orfs[id], gp = GPOS[id], aa = o.aa || "", nt = o.nt || "";
        var ncod = Math.floor(gp.length / 3);
        for (var k = 0; k < ncod; k++) {
            var p0 = gp[3 * k], p1 = gp[3 * k + 1], p2 = gp[3 * k + 2];
            if (Math.max(p0, p1, p2) < viewG[0] - 2 || Math.min(p0, p1, p2) > viewG[1] + 2) continue;
            var res = k < aa.length ? aa.charAt(k) : "*";
            var isStart = k === 0, isStop = k >= aa.length;
            var fill = isStart ? "#1f6b45" : isStop ? "#a5322c" : "#5b6b82";
            var bg = isStart ? "#d9efe1" : isStop ? "#f3dcda" : "#eef2f7";
            var inView = [ p0, p1, p2 ].filter(function(pp) {
                return pp >= viewG[0] - 2 && pp <= viewG[1] + 2;
            }).map(function(pp) {
                return x(pp + .5);
            });
            var cx = inView.length ? inView.reduce(function(s, v) {
                return s + v;
            }, 0) / inView.length : x(p1 + .5);
            g.append("line").attr("x1", x(Math.min(p0, p2) + 0)).attr("x2", x(Math.max(p0, p2) + 1)).attr("y1", y + 11).attr("y2", y + 11).attr("stroke", fill).attr("stroke-width", 1).attr("opacity", .5);
            g.append("rect").attr("x", cx - 5).attr("y", y).attr("width", 10).attr("height", 12).attr("rx", 2).attr("fill", bg);
            g.append("text").attr("x", cx).attr("y", y + 9.5).attr("text-anchor", "middle").attr("font-size", 9.5).attr("font-family", MONO).attr("font-weight", 600).attr("fill", fill).text(res);
        }
        hdr.append("text").attr("x", 12).attr("y", y + 10).attr("font-size", 9).attr("font-weight", 600).attr("fill", "#5f6a66").text("ORF" + id);
    }
    function isCarrierOf(orfId, t) {
        return !!orfId && META.orfs[orfId].carriers.indexOf(t) >= 0;
    }
    function rowOpacityFor(t) {
        return selectedOrf && !isCarrierOf(selectedOrf, t) ? .35 : 1;
    }
    function clipped(g) {
        return g.append("g").attr("clip-path", "url(#plot-clip)");
    }
    function renderStack() {
        svg.selectAll("*").remove();
        svg.append("defs").append("clipPath").attr("id", "plot-clip").append("rect").attr("x", plotL()).attr("y", 0).attr("width", Math.max(1, plotR() - plotL())).attr("height", H);
        gMain = svg.append("g");
        gMain.append("rect").attr("class", "deselect-catcher").attr("x", 0).attr("y", 0).attr("width", W).attr("height", H).attr("fill", "transparent").on("click", function() {
            pickOrf(null);
        });
        gMain.append("rect").attr("class", "gutter-bg").attr("x", 0).attr("y", 0).attr("width", plotL() - 1).attr("height", H).attr("fill", "#ffffff").style("pointer-events", "none");
        var gGrid = clipped(gMain).style("pointer-events", "none");
        x.ticks(8).forEach(function(v) {
            var xp = x(v);
            if (xp < plotL() || xp > plotR()) return;
            gGrid.append("line").attr("x1", xp).attr("x2", xp).attr("y1", 0).attr("y2", H).attr("stroke", "#eef1f3");
        });
        if (selectedOrf) {
            var sp = META.orfs[selectedOrf].span;
            gGrid.append("rect").attr("x", x(sp[0])).attr("y", 0).attr("width", Math.max(1.5, x(sp[1] + 1) - x(sp[0]))).attr("height", H).attr("fill", PAL.nmyc).attr("opacity", .06);
        }
        gMain.append("line").attr("x1", plotL() - 6).attr("x2", plotL() - 6).attr("y1", 0).attr("y2", H).attr("stroke", "#e3e7e9").style("pointer-events", "none");
        ROWS.forEach(function(r) {
            if (r.kind === "tx") drawTx(r); else if (r.kind === "grp") drawGroup(r); else if (r.kind === "orf") drawOrf(r); else if (r.kind === "region") drawRegion(r); else if (r.kind === "shared") drawShared(r);
        });
        drawMsGroups();
    }
    function drawMsGroups() {}
    function msGroupTip() {
        return "Experimental support has not been assigned in this reference build.";
    }
    function dmycnTx() {
        return TX_ORDER.find(function(t) {
            var p = META.transcripts[t].protein_sequence;
            return p && p !== META.transcripts["201"].protein_sequence;
        }) || null;
    }
    function segFill(t, type) {
        if (type === "cds") return t === dmycnTx() ? PAL.dmycn : PAL.nmyc;
        if (type === "ri") return PAL.ri_line;
        return PAL.utr;
    }
    function drawOrfHighlight(gClip, orfId, mid) {
        var o = META.orfs[orfId];
        var b0 = x(o.blocks[0][0]), b1 = x(o.blocks[o.blocks.length - 1][1] + 1);
        gClip.append("line").attr("x1", b0).attr("x2", b1).attr("y1", mid).attr("y2", mid).attr("stroke", PAL.nmyc).attr("stroke-width", 1).attr("opacity", .55).style("pointer-events", "none");
        o.blocks.forEach(function(b) {
            var xa = x(b[0]), xb = x(b[1] + 1), w = Math.max(2, xb - xa);
            gClip.append("rect").attr("class", "orf-hl-halo").attr("x", xa - 1.5).attr("y", mid - 8).attr("width", w + 3).attr("height", 16).attr("rx", 3).attr("fill", o.color).attr("opacity", .3);
            gClip.append("rect").attr("class", "orf-hl").attr("x", xa).attr("y", mid - 6).attr("width", w).attr("height", 12).attr("rx", 2.5).attr("fill", o.color).attr("stroke", PAL.nmyc).attr("stroke-width", 1.6).on("click", function(ev) {
                ev.stopPropagation();
                pickOrf(orfId);
            }).on("mouseenter", function(ev) {
                tip(ev, orfTip(orfId));
            }).on("mousemove", moveTip).on("mouseleave", hideTip);
        });
    }
    function scrollCarrierIntoView(orfId) {
        var sc = $("stack-scroll");
        if (!sc || !orfId) return;
        var cs = META.orfs[orfId].carriers;
        var row = null;
        for (var i = 0; i < ROWS.length; i++) {
            if (ROWS[i].kind === "tx" && cs.indexOf(ROWS[i].t) >= 0) {
                row = ROWS[i];
                break;
            }
        }
        if (!row) return;
        var top = row.y, bottom = row.y + TXH, vt = sc.scrollTop, vh = sc.clientHeight || 400;
        if (top < vt) sc.scrollTop = Math.max(0, top - 8); else if (bottom > vt + vh) sc.scrollTop = bottom - vh + 8;
    }
    function drawTx(r) {
        var t = r.t, tx = META.transcripts[t], mid = r.y + TXH / 2, op = rowOpacityFor(t);
        var g = gMain.append("g").attr("class", "tx-row").attr("data-tx", t).attr("opacity", op);
        var gc = clipped(g);
        g.append("rect").attr("class", "gutter-hit").attr("x", 0).attr("y", r.y).attr("width", GUTTER - 8).attr("height", TXH).attr("fill", "transparent").on("click", function() {
            pickTx(t);
        }).on("mouseenter", function(ev) {
            tip(ev, "<b>" + esc(LAB.label(t)) + "</b><br>Laboratory alias: " + esc(LAB.alias(t)) + "<br>Reference: " + esc(txRel(t)) + "<br><span class='mono'>" + esc(txAcc(t)) + "</span><br><span class='dim'>" + esc(txSrc(t)) + "<br>" + orfIdsOfTx(t).length + " ORFs<br>click to " + (state[t].open ? "collapse" : "expand") + "</span>");
        }).on("mousemove", moveTip).on("mouseleave", hideTip);
        var sh = shortLabels(), nOrf = orfIdsOfTx(t).length, coding = txCoding(t);
        var yy = r.y + 16;
        txLines(t).forEach(function(L) {
            var tn = g.append("text").attr("class", L.cls).attr("x", 14).attr("y", yy).attr("font-size", L.fs).attr("fill", L.cls === "tx-label" ? selectedTx === t ? PAL.nmyc : "#222a2e" : L.cls === "tx-acc" ? txOffV49(t) ? "#a5322c" : "#5f6a66" : L.cls === "tx-class" ? coding === "non-coding" ? "#a5322c" : "#5f6a66" : "#5f6a66");
            if (L.fw) tn.attr("font-weight", L.fw);
            if (L.mono) tn.attr("font-family", MONO);
            tn.node().textContent = L.txt;
            yy += L.fs + 2.5;
        });
        var zx = GUTTER - 16;
        g.append("text").attr("class", "tx-chev").attr("x", zx).attr("y", r.y + 20).attr("font-size", 12).attr("fill", "#5f6a66").attr("text-anchor", "middle").text(state[t].open ? "▾" : "▸");
        var x0e = x(tx.exons[0][0]), x1e = x(tx.exons[tx.exons.length - 1][1] + 1);
        gc.append("line").attr("class", "backbone").attr("x1", x0e).attr("x2", x1e).attr("y1", mid).attr("y2", mid).attr("stroke", BACKBONE).attr("stroke-width", 1.2);
        chevrons(gc, x0e, x1e, mid, BACKBONE);
        tx.exons.forEach(function(exon, ei) {
            exonSegments(t, exon).forEach(function(s) {
                var h = s.type === "cds" ? THICK : s.type === "ri" ? RIH : THIN, xa = x(s.a), xb = x(s.b + 1), w = Math.max(1.2, xb - xa);
                gc.append("rect").attr("class", "exon").attr("x", xa).attr("y", mid - h / 2).attr("width", w).attr("height", h).attr("rx", 1.5).attr("fill", segFill(t, s.type)).attr("role", "button").attr("tabindex", 0).attr("aria-label", LAB.accessible(t) + ", exon " + (ei + 1) + ", chr2:" + exon[0] + "–" + exon[1] + ", " + s.type).on("keydown", function(ev) {
                    if (ev.key === "Enter" || ev.key === " ") {
                        ev.preventDefault();
                        showExon(t, exon, ei);
                        $("details-close").focus();
                    }
                }).on("click", function(ev) {
                    ev.stopPropagation();
                    showExon(t, exon, ei);
                }).on("mouseenter", function(ev) {
                    tip(ev, exonTip(t, exon, ei, s));
                }).on("mousemove", moveTip).on("mouseleave", hideTip);
            });
        });
        if (isCarrierOf(selectedOrf, t)) drawOrfHighlight(gc, selectedOrf, mid);
    }
    function drawGroup(r) {
        var t = gMain.append("text").attr("class", "grp-lbl").attr("x", 20).attr("y", r.y + 14).attr("font-size", 10).attr("font-weight", 600).attr("letter-spacing", ".5px").attr("fill", r.uniq ? PAL.nmyc : "#5f6a66").text(r.label.toUpperCase());
        if (r.uniq) t.append("title").text("within the "+TX_ORDER.length+" transcripts shown");
    }
    function drawShared(r) {
        var t = r.t, g = gMain.append("g");
        g.append("rect").attr("class", "shared-toggle-hit").attr("x", 20).attr("y", r.y + 3).attr("width", GUTTER - 40).attr("height", SHH - 8).attr("rx", 5).attr("fill", "#f2f4f6").attr("stroke", "#d5dbe2").attr("stroke-dasharray", "3 3").attr("role", "button").attr("tabindex", 0).attr("aria-label", LAB.accessible(t) + ", toggle shared ORFs").attr("aria-expanded", String(state[t].sharedOpen)).on("keydown", function(ev) {
            if (ev.key === "Enter" || ev.key === " ") {
                ev.preventDefault();
                state[t].sharedOpen = !state[t].sharedOpen;
                rebuild();
                var n = document.querySelector('.gutter-hit[aria-expanded="true"]');
                if (n) n.focus();
            }
        }).on("click", function() {
            state[t].sharedOpen = !state[t].sharedOpen;
            rebuild();
        });
        g.append("text").attr("class", "shared-toggle-lbl").attr("x", 25).attr("y", r.y + 15).attr("font-size", 10).attr("font-weight", 600).attr("fill", PAL.nmyc).text("Shared ORFs").append("tspan").attr("x", 25).attr("dy", 14).text((state[t].sharedOpen ? "▾ Hide " : "▸ Show ") + r.n);
    }
    function regionOf() {
        return null;
    }
    function regionRec(rid) {
        return (META.regions || {})[rid] || null;
    }
    function regionCollapses() {
        return false;
    }
    function regionMembers(rid) {
        var R = regionRec(rid);
        return R ? R.members.map(String) : [];
    }
    function eventsOf(ids) {
        return ids.map(function(id) {
            return {
                kind: "orf",
                id: id
            };
        });
    }
    function regionFrames(rid) {
        var R = regionRec(rid);
        if (!R) return 1;
        var f = {};
        R.members.forEach(function(i) {
            var o = META.orfs[String(i)];
            if (o) f[o.span[1] % 3] = 1;
        });
        return Object.keys(f).length;
    }
    function landingOrfs() {
        return Object.keys(META.orfs).filter(function(id) {
            return META.orfs[id].annotation_class === "ANNOTATED_CDS";
        });
    }
    function renderLanding() {
        return "<h2>Current reference</h2><p>" + esc(META.reference_build) + "</p><p>Complete GENCODE MYCN set plus the reviewed RefSeq MYCN-XM addition. ORFs are independently enumerated from current sequences; local ORF numbers do not correspond to the historical viewer.</p><p>Annotations and predictions do not establish expression or function. See Reference and evidence for source versions, validation, and downloadable data.</p>";
    }
    function drawRegion() {}
    function regionTip() {
        return "No experimental region assignments in the current build.";
    }
    var _mrule = null;
    function textW(txt, fs, fw, ff) {
        if (!_mrule || !_mrule.node() || !_mrule.node().isConnected) {
            _mrule = svg.append("text").attr("class", "measure-only").attr("x", -9999).attr("y", -9999).attr("visibility", "hidden");
        }
        _mrule.attr("font-size", fs).attr("font-weight", fw || null).attr("font-family", ff || FONT);
        _mrule.node().textContent = txt;
        var w = _mrule.node().getComputedTextLength ? _mrule.node().getComputedTextLength() : 0;
        if (w === 0 && txt && txt.length) w = txt.length * fs * .62;
        return w;
    }
    function packTokens(tokens, avail, fs, fw, ff) {
        var lines = [], cur = "";
        tokens.forEach(function(tk) {
            var cand = cur ? cur + " · " + tk : tk;
            if (textW(cand, fs, fw, ff) <= avail) cur = cand; else {
                if (cur) lines.push(cur);
                cur = tk;
            }
        });
        if (cur) lines.push(cur);
        return lines.length ? lines : [ "" ];
    }
    function txLines(t) {
        var avail = GUTTER - 22 - 14;
        var coding = txCoding(t), nOrf = orfIdsOfTx(t).length, nEv = eventCountOfTx(t);
        var out = [];
        out.push({
            txt: LAB.canonical(t),
            fs: 13,
            fw: 600,
            cls: "tx-label"
        });
        if (LAB.alias(t)) out.push({
            txt: "Lab alias: " + LAB.alias(t),
            fs: 11,
            cls: "tx-alias"
        });
        packTokens([ coding === "non-coding" ? "no annotated CDS" : "annotated CDS" ], avail, 10.5).forEach(function(l) {
            out.push({
                txt: l,
                fs: 10.5,
                cls: "tx-class"
            });
        });
        packTokens((txAcc(t) + txTagShort(t)).split(" · "), avail, 9.5, null, MONO).forEach(function(l) {
            out.push({
                txt: l,
                fs: 9.5,
                cls: "tx-acc",
                mono: true
            });
        });
        packTokens([ nOrf + " ORFs" ], avail, 10.5).forEach(function(l) {
            out.push({
                txt: l,
                fs: 10.5,
                cls: "tx-meta"
            });
        });
        return out;
    }
    function orfLines(id) {
        var o = META.orfs[id], avail = GUTTER - 8 - 42, nm = NAME[id] ? " " + NAME[id] : "";
        var out = [];
        packTokens([ "ORF" + id + nm ], avail, 11, 600).forEach(function(l) {
            out.push({
                txt: l,
                fs: 11,
                fw: 600,
                cls: "orf-lbl"
            });
        });
        packTokens([ o.aa_len + " aa", carrierCount(id) + "c" ], avail, 10).forEach(function(l) {
            out.push({
                txt: l,
                fs: 10,
                cls: "orf-tag"
            });
        });
        return out;
    }
    function orfRowH() {
        var mx = 1;
        Object.keys(META.orfs).forEach(function(id) {
            mx = Math.max(mx, orfLines(id).length);
        });
        return 6 + mx * 12 + 6;
    }
    function txRowH() {
        var mx = 1;
        TX_ORDER.forEach(function(t) {
            mx = Math.max(mx, txLines(t).length);
        });
        return 12 + mx * 13 + 6;
    }
    function drawOrf(r) {
        var id = r.id, o = META.orfs[id], mid = r.y + ORFH / 2;
        var g = gMain.append("g").attr("class", "orf-g").attr("data-id", id);
        var gc = clipped(g);
        g.append("rect").attr("class", "orf-hit").attr("x", 22).attr("y", r.y).attr("width", GUTTER - 30).attr("height", ORFH).attr("fill", "transparent").style("cursor", "pointer").on("click", function() {
            pickOrf(id);
        }).on("mouseenter", function(ev) {
            tip(ev, orfTip(id));
        }).on("mousemove", moveTip).on("mouseleave", hideTip);
        g.append("rect").attr("class", "orf-sw").attr("x", 24).attr("y", r.y + 5).attr("width", 10).attr("height", 10).attr("rx", 2).attr("fill", o.color).attr("stroke", msDetected(o) ? PAL.evid : "none").attr("stroke-width", msDetected(o) ? 1.6 : 0);
        var yy2 = r.y + 13;
        orfLines(id).forEach(function(L) {
            var tn2 = g.append("text").attr("class", L.cls).attr("x", 42).attr("y", yy2).attr("font-size", L.fs).attr("fill", L.cls === "orf-lbl" ? "#222a2e" : "#5f6a66").style("cursor", "pointer").on("click", function() {
                pickOrf(id);
            });
            if (L.fw) tn2.attr("font-weight", L.fw);
            tn2.node().textContent = L.txt;
            yy2 += L.fs + 2;
        });
        var b0 = x(o.blocks[0][0]), b1 = x(o.blocks[o.blocks.length - 1][1] + 1);
        gc.append("line").attr("x1", b0).attr("x2", b1).attr("y1", mid).attr("y2", mid).attr("stroke", "#c2ccc9").attr("stroke-width", 1);
        chevrons(gc, b0, b1, mid, "#c2ccc9");
        o.blocks.forEach(function(b) {
            var xa = x(b[0]), xb = x(b[1] + 1), w = Math.max(1.5, xb - xa);
            gc.append("rect").attr("class", "orf-bar").attr("x", xa).attr("y", mid - ORF_THICK / 2).attr("width", w).attr("height", ORF_THICK).attr("rx", 2).attr("fill", o.color).on("click", function() {
                pickOrf(id);
            }).on("mouseenter", function(ev) {
                tip(ev, orfTip(id));
            }).on("mousemove", moveTip).on("mouseleave", hideTip);
        });
        if (selectedOrf === id) {
            gc.append("rect").attr("class", "orf-sel").attr("x", b0 - 3).attr("y", mid - ORF_THICK / 2 - 3).attr("width", b1 - b0 + 6).attr("height", ORF_THICK + 6).attr("rx", 4).attr("fill", "none").attr("stroke", PAL.nmyc).attr("stroke-width", 1.6);
        }
    }
    function setNow(name, g) {
        selNowName = name;
        drawNow(name, g || viewG);
    }
    var selNowName = "MYCN locus";
    function updateNow() {
        drawNow(selNowName, viewG);
    }
    function drawNow(name, g) {
        $("ns-name").textContent = name;
        $("ns-coord").textContent = "chr2:" + fmt(g[0]) + "-" + fmt(g[1]);
    }
    function paddedSpan(sp) {
        var pad = Math.max(30, (sp[1] - sp[0]) * .12);
        return [ clamp(sp[0] - pad, DMIN, DMAX), clamp(sp[1] + pad, DMIN, DMAX) ];
    }
    function pickTx(t) {
        if (selectedOrf != null) pickOrf(null);
        state[t].open = !state[t].open;
        if (!state[t].open) state[t].sharedOpen = false;
        selectedTx = state[t].open ? t : null;
        rebuild();
        setNow(state[t].open ? LAB.label(t) : "MYCN locus", null);
        renderFeatureList();
    }
    function pickRegion(rid) {
        rid = rid == null ? null : String(rid);
        if (rid != null && rid === selectedRegion) rid = null;
        selectedRegion = rid;
        selectedOrf = null;
        if (rid == null) {
            railHasContent = false;
            $("rail-body").innerHTML = renderLanding();
            wireCard();
        } else {
            $("rail-body").innerHTML = renderRegionCard(rid);
            wireCard();
            railHasContent = true;
            openRight(true);
        }
        rebuild();
    }
    function renderRegionCard() {
        return renderLanding();
    }
    function pickOrf(id) {
        id = id == null ? null : String(id);
        if (id != null && id === selectedOrf) id = null;
        selectedOrf = id;
        if (id != null) selectedRegion = null;
        var target = viewG.slice();
        if (id == null) {
            railHasContent = false;
            $("rail-body").innerHTML = renderLanding();
            wireCard();
            openRight(false);
            setNow(selectedTx ? LAB.label(selectedTx) : "MYCN locus", selectedTx ? txSpan(selectedTx) : null);
        } else {
            var o = META.orfs[id], sp = o.span;
            $("rail-body").innerHTML = renderCard(id, o);
            wireCard();
            railHasContent = true;
            var fc = TX_ORDER.filter(function(tt) {
                return carriers(id).indexOf(tt) >= 0;
            })[0];
            if (fc) state[fc].open = true;
            var visible = sp[1] >= viewG[0] && sp[0] <= viewG[1];
            if (!visible) target = paddedSpan(sp);
            setNow("ORF" + id + (NAME[id] ? " · " + NAME[id] : ""), visible ? null : sp);
            openRight(true);
        }
        rebuild();
        applyRange(target[0], target[1], false);
        scrollCarrierIntoView(id);
        markFeatureActive();
    }
    function showExon(t, exon, ei) {
        var cds = META.transcripts[t].cds || [], isCds = false, phase = null;
        cds.forEach(function(c, ci) {
            if (!(c[1] < exon[0] || c[0] > exon[1])) {
                isCds = true;
                if (phase === null) phase = CDS_PHASE[t][ci];
            }
        });
        var len = exon[1] - exon[0] + 1, riHit = (RETAINED[t] || []).some(function(iv) {
            return !(iv[1] < exon[0] || iv[0] > exon[1]);
        });
        var h = '<div class="exon-detail"><h2>Exon detail</h2>';
        h += '<p class="exon-title">' + esc(LAB.label(t)) + " · exon " + (ei + 1) + '</p><table class="props">';
        h += '<tr><td class="k">Coordinates</td><td class="v">chr2:' + fmt(exon[0]) + "–" + fmt(exon[1]) + "</td></tr>";
        h += '<tr><td class="k">Length</td><td class="v">' + fmt(len) + " bp</td></tr>";
        h += '<tr><td class="k">Content</td><td class="v">' + (isCds ? "contains coding sequence (may include UTR)" : txCoding(t) === "non-coding" ? "non-coding exon" : "untranslated exon") + "</td></tr>";
        if (isCds) h += '<tr><td class="k">CDS phase (+ strand)</td><td class="v">' + phase + "</td></tr>";
        if (riHit) h += '<tr><td class="k">Retained intron</td><td class="v">yes (' + esc(t) + ")</td></tr>";
        h += '<tr><td class="k">Strand</td><td class="v">' + esc(META.transcripts[t].strand) + "</td></tr>";
        var _x = META.transcripts[t] || {};
        h += propRow("Laboratory alias", LAB.alias(t));
        h += transcriptProperties(t);
        h += '<tr><td class="k">Transcript ID</td><td class="v mono">' + esc(txAcc(t)) + "</td></tr>";
        h += '<tr><td class="k">Source</td><td class="v">' + esc(txSrc(t)) + "</td></tr>";
        h += '<tr><td class="k">Source release</td><td class="v">' + esc(txRel(t)) + "</td></tr>";
        if (_x.refseq_variant) h += '<tr><td class="k">Variant</td><td class="v">' + esc(_x.refseq_variant) + "</td></tr>";
        if (_x.refseq_cds) h += '<tr><td class="k">CDS (mRNA)</td><td class="v mono">' + esc(_x.refseq_cds) + "</td></tr>";
        if (_x.refseq_protein) h += '<tr><td class="k">Protein</td><td class="v mono">' + esc(_x.refseq_protein) + "</td></tr>";
        if (_x.refseq_length_bp) h += '<tr><td class="k">Length</td><td class="v mono">' + _x.refseq_length_bp + " bp · " + _x.refseq_exons + " exons</td></tr>";
        if (_x.refseq_obsolete_id) h += '<tr><td class="k">Superseded model</td><td class="v mono">' + esc(_x.refseq_obsolete_id) + ' <span class="obs">OBSOLETE</span></td></tr>';
        if (_x.refseq_absorbed && _x.refseq_absorbed.length) h += '<tr><td class="k">Absorbed by the NM</td><td class="v mono">' + esc(_x.refseq_absorbed.join(" · ")) + "</td></tr>";
        var _m = _x.mane_select;
        if (_m && _m.is_mane_select) {
            h += '<tr><td class="k">MANE Select</td><td class="v mono">' + esc(_m.ensembl_transcript) + " · " + esc(_m.refseq_transcript) + "</td></tr>";
            h += '<tr><td class="k">MANE proteins</td><td class="v mono">' + esc(_m.refseq_protein) + " · " + esc(_m.ensembl_protein) + "</td></tr>";
        }
        h += "</table>";
        if (_m && _m.is_mane_select) {
            h += '<p class="txnote mane"><b>MANE Select.</b> ' + esc(_m.note) + "</p>";
            h += '<p class="txnote warn">' + esc(_m.version_discrepancy) + "</p>";
        }
        if (_x.model_verified) h += '<p class="txnote mane">' + esc(_x.model_verified) + "</p>";
        if (_x.refseq_current_note) h += '<p class="txnote">' + esc(_x.refseq_current_note) + "</p>";
        if (_x.refseq_obsolete_note) h += '<p class="txnote warn">' + esc(_x.refseq_obsolete_note) + "</p>";
        if (_x.id_note) h += '<p class="txnote warn">' + esc(_x.id_note) + "</p>";
        h += '<p class="rail-empty">Click an ORF (expand the transcript) to see its protein detail.</p></div>';
        $("rail-body").innerHTML = h;
        railHasContent = true;
        openRight(true);
        setNow(LAB.label(t) + " exon " + (ei + 1), exon);
        var target = viewG.slice();
        rebuild();
        applyRange(target[0], target[1], false);
    }
    function tip(ev, html) {
        var el = $("tip");
        el.innerHTML = html;
        el.hidden = false;
        moveTip(ev);
    }
    function moveTip(ev) {
        var el = $("tip");
        if (el.hidden) return;
        var vw = window.innerWidth || 1200, vh = window.innerHeight || 800, tw = el.offsetWidth || 180, th = el.offsetHeight || 40;
        var lx = ev.clientX + 14, ty = ev.clientY + 14;
        if (lx + tw + 8 > vw) lx = Math.max(8, ev.clientX - tw - 14);
        if (ty + th + 8 > vh) ty = Math.max(8, ev.clientY - th - 14);
        el.style.left = lx + "px";
        el.style.top = ty + "px";
    }
    function hideTip() {
        $("tip").hidden = true;
    }
    function exonTip(t, exon, ei, s) {
        var kind = s.type === "cds" ? "CDS" : s.type === "ri" ? "retained intron" : txCoding(t) === "non-coding" ? "no annotated CDS" : "UTR";
        return "<b>" + esc(LAB.label(t)) + " exon " + (ei + 1) + "</b> · " + kind + "<br>Laboratory alias: " + esc(LAB.alias(t)) + "<br>" + esc(txAcc(t)) + "<br>" + esc(txSrc(t)) + "<br>Reference: " + esc(txRel(t)) + "<br><span class='dim'>chr2:" + fmt(exon[0]) + "-" + fmt(exon[1]) + " · " + fmt(exon[1] - exon[0] + 1) + " bp</span>";
    }
    function orfTip(id) {
        var o = META.orfs[id], nm = NAME[id] ? " · " + NAME[id] : "";
        return "<b>ORF" + id + nm + "</b> · " + esc(o.cls) + "<br><span class='dim'>" + o.aa_len + " aa · " + carrierCount(id) + " carrier" + (carrierCount(id) === 1 ? "" : "s") + " · click for detail</span>";
    }
    function setLens(l) {
        if (l === lens) return;
        lens = l;
        $("lens-tx").classList.toggle("active", l === "tx");
        $("lens-orf").classList.toggle("active", l === "orf");
        $("cue").firstChild && ($("cue").childNodes[0].nodeValue = l === "tx" ? "Select a transcript to explore its annotated and predicted ORFs." : "Click an ORF to view it on the stack.");
        renderFeatureList();
        var sq = $("feature-search").value;
        if (sq) filterList(sq);
    }
    function renderFeatureList() {
        var active = document.activeElement, key = active && active.getAttribute("data-tx");
        if (lens === "orf") renderOrfList(); else renderTxList();
        markFeatureActive();
        if (key) {
            var next = document.querySelector('.trow[data-tx="' + key + '"]');
            if (next) next.focus();
        }
        if ($("feature-search").value) filterList($("feature-search").value);
    }
    function renderTxList() {
        novelOpen = false;
        var host = $("feature-list"), h = '<div class="fgroup">';
        TX_ORDER.forEach(function(t) {
            var n = orfIdsOfTx(t).length, coding = txCoding(t), nEv = eventCountOfTx(t);
            h += '<button class="trow" data-tx="' + t + '" title="' + esc(LAB.label(t)) + '" aria-label="' + esc(LAB.accessible(t) + ", " + txAcc(t)) + '" aria-expanded="' + state[t].open + '">' + '<span class="tchev">' + (state[t].open ? "▾" : "▸") + "</span>" + '<span class="tname">' + esc(LAB.label(t)) + "</span>" + '<span class="tacc' + (txOffV49(t) ? " toff" : "") + '">' + esc(txAcc(t)) + esc(txTag(t)) + "</span>" + '<span class="tmeta">' + n + ' ORFs · <span class="' + (coding === "non-coding" ? "tnc" : "") + '">' + (coding === "non-coding" ? "no annotated CDS" : "annotated CDS") + "</span></span></button>";
        });
        h += "</div>";
        host.innerHTML = h;
        Array.prototype.forEach.call(host.querySelectorAll(".trow"), function(btn) {
            btn.addEventListener("click", function() {
                pickTx(btn.getAttribute("data-tx"));
            });
        });
    }
    function orfRowHTML(id) {
        var o = META.orfs[id];
        return '<button class="frow" data-id="' + id + '"><span class="sw" style="background:' + o.color + '"></span><span class="ftxt">ORF' + id + " · " + o.aa_len + " aa · " + esc(o.annotation_class) + '</span><span class="ntag">' + carrierCount(id) + "c</span></button>";
    }
    function renderOrfList() {
        novelOpen = false;
        var host = $("feature-list"), h = "";
        GROUPS.forEach(function(gr) {
            var ids = orfsOfClass(gr.key);
            if (!ids.length) return;
            h += '<div class="fgroup"><div class="fghead">' + esc(gr.title) + "</div>" + ids.map(orfRowHTML).join("") + "</div>";
        });
        var novelIds = NOVEL.reduce(function(a, n) {
            return a.concat(orfsOfClass(n.key));
        }, []);
        h += '<div class="fgroup novelwrap"><button id="novel-toggle" class="novel-toggle" aria-expanded="false">' + '<span class="chev">▸</span> Show ' + novelIds.length + ' predicted ORFs</button><div id="novel-body" class="novel-body" hidden>';
        NOVEL.forEach(function(n) {
            var ids = orfsOfClass(n.key);
            h += '<div class="fsub"><div class="fghead sub">' + esc(n.title) + " (" + ids.length + ")</div>" + ids.map(orfRowHTML).join("") + "</div>";
        });
        h += "</div></div>";
        host.innerHTML = h;
        Array.prototype.forEach.call(host.querySelectorAll(".frow"), function(btn) {
            btn.addEventListener("click", function() {
                pickOrf(btn.getAttribute("data-id"));
            });
        });
        var ntog = $("novel-toggle");
        if (ntog) ntog.addEventListener("click", function() {
            setNovel(!novelOpen);
        });
    }
    function setNovel(open) {
        novelOpen = open;
        var b = $("novel-body");
        if (!b) return;
        b.hidden = !open;
        var t = $("novel-toggle");
        t.setAttribute("aria-expanded", open ? "true" : "false");
        t.querySelector(".chev").textContent = open ? "▾" : "▸";
    }
    function filterList(q) {
        q = (q || "").trim().toLowerCase();
        var host = $("feature-list"), anyNovel = false;
        Array.prototype.forEach.call(host.querySelectorAll(".frow"), function(btn) {
            var id = btn.getAttribute("data-id"), o = META.orfs[id];
            var hay = [ "ORF" + id, NAME[id] || "", o.cls, (o.carriers || []).join(" "), o.evidence_class, o.annotation_class, o.acc || "" ].join(" ").toLowerCase();
            var hit = /^(?:tv\d+\+?|nc)$/.test(q) ? (o.carriers || []).some(function(t) {
                return LAB.matches(t, txAcc(t), q);
            }) : !q || hay.indexOf(q) >= 0;
            btn.hidden = !hit;
            if (hit && q && btn.closest(".novel-body")) anyNovel = true;
        });
        Array.prototype.forEach.call(host.querySelectorAll(".trow"), function(btn) {
            var t = btn.getAttribute("data-tx");
            btn.hidden = !LAB.matches(t, txAcc(t), q, [ txCoding(t), META.transcripts[t].protein_accession, META.transcripts[t].description, META.transcripts[t].biotype ].join(" "));
        });
        if (q && anyNovel && !novelOpen) setNovel(true);
        var visible = Array.prototype.filter.call(host.querySelectorAll(".trow,.frow"), function(b) {
            return !b.hidden;
        }).length;
        $("search-status").textContent = visible ? visible + " matching features" : "No matching features";
    }
    function markFeatureActive() {
        Array.prototype.forEach.call(document.querySelectorAll(".frow.active,.trow.active"), function(b) {
            b.classList.remove("active");
        });
        if (lens === "orf" && selectedOrf) {
            var b = document.querySelector('.frow[data-id="' + selectedOrf + '"]');
            if (b) b.classList.add("active");
        }
        if (lens === "tx" && selectedTx) {
            var t = document.querySelector('.trow[data-tx="' + selectedTx + '"]');
            if (t) t.classList.add("active");
        }
    }
    function wireCard() {
        var body = $("rail-body");
        Array.prototype.forEach.call(body.querySelectorAll("button.copy"), function(btn) {
            btn.addEventListener("click", function() {
                copySeq(btn);
            });
        });
        Array.prototype.forEach.call(body.querySelectorAll("button.memb"), function(btn) {
            btn.addEventListener("click", function() {
                pickOrf(btn.getAttribute("data-orf"));
            });
        });
        Array.prototype.forEach.call(body.querySelectorAll("button.tile"), function(btn) {
            btn.addEventListener("click", function() {
                pickOrf(btn.getAttribute("data-orf"));
            });
        });
    }
    function copySeq(btn) {
        var txt = btn.getAttribute("data-seq") || "", old = btn.getAttribute("data-label") || btn.textContent;
        btn.setAttribute("data-label", old);
        function ok() {
            btn.textContent = "copied";
            btn.classList.add("ok");
            setTimeout(function() {
                btn.textContent = old;
                btn.classList.remove("ok");
            }, 1100);
        }
        function fail() {
            btn.textContent = "copy failed";
            btn.classList.add("err");
            setTimeout(function() {
                btn.textContent = old;
                btn.classList.remove("err");
            }, 1400);
        }
        if (navigator.clipboard && navigator.clipboard.writeText) {
            navigator.clipboard.writeText(txt).then(ok, function() {
                if (!legacyCopy(txt)) fail(); else ok();
            });
        } else {
            if (legacyCopy(txt)) ok(); else fail();
        }
    }
    function legacyCopy(txt) {
        try {
            var ta = document.createElement("textarea");
            ta.value = txt;
            ta.style.position = "fixed";
            ta.style.opacity = "0";
            document.body.appendChild(ta);
            ta.select();
            var okc = document.execCommand && document.execCommand("copy");
            document.body.removeChild(ta);
            return !!okc;
        } catch (e) {
            return false;
        }
    }
    function badgesFor(id, o) {
        return '<span class="badge">' + esc(o.evidence_class) + "</span>";
    }
    function propRow(k, v) {
        return v === null || v === undefined || v === "" ? "" : '<tr><td class="k">' + esc(k) + '</td><td class="v">' + esc(v) + "</td></tr>";
    }
    function propRowNA(k, reason) {
        return '<tr><td class="k">' + esc(k) + '</td><td class="v na">' + esc(reason) + "</td></tr>";
    }
    function slug(s) {
        return String(s).toLowerCase().replace(/[^a-z]+/g, "-").replace(/^-|-$/g, "");
    }
    function evidenceBlock(o) {
        return "<p>" + esc(o.source) + ". No sample-specific expression or function is inferred.</p>";
    }
    function codonStrip(nt, aa) {
        nt = String(nt || "");
        aa = String(aa || "");
        var out = [];
        for (var i = 0; i < nt.length; i += 3) {
            var codon = nt.substr(i, 3), idx = i / 3, res = idx < aa.length ? aa[idx] : "*";
            var cls = i === 0 ? "start" : idx >= aa.length ? "stop" : "";
            out.push('<span class="codon ' + cls + '"><b>' + esc(codon) + "</b><i>" + esc(res) + "</i></span>");
        }
        return out.join("");
    }
    function msSummary() {
        return "Not assigned: no current record-level experimental evidence join.";
    }
    function msCaveat() {
        return "";
    }
    function renderCard(id, o) {
        var h = "<h2>ORF " + esc(id) + "</h2><p>" + esc(o.build_orf_id) + "</p><p>" + esc(o.annotation_class) + " · " + esc(o.evidence_class) + "</p>";
        h += "<p>" + esc(o.source) + '</p><p>ORF identity is specific to this reference build. Annotation and prediction do not establish expression or function.</p><table class="props">';
        [ [ "Protein accession", o.acc || "No annotated protein accession" ], [ "Length", o.aa_len + " aa" ], [ "Genomic blocks (inclusive)", o.blocks.map(function(b) {
            return b.join("–");
        }).join("; ") ], [ "Carrier transcripts", o.carriers.map(LAB.label).join(", ") ], [ "Relationships", o.location ], [ "Start / stop", o.start_codon + " / " + o.stop_codon ], [ "Molecular weight (Da)", o.mw ], [ "Theoretical pI", o.pI ], [ "GRAVY", o.gravy ], [ "Instability index (prediction)", o.instab ], [ "Net charge (pH 7)", o.netq ], [ "Property algorithm", "Biopython " + o.property_provenance.version + " ProteinAnalysis" ], [ "Input protein SHA-256", o.protein_sequence_sha256 ] ].forEach(function(row) {
            h += propRow(row[0], row[1]);
        });
        h += "</table>";
        h += "<h3>Transcript-coordinate mappings</h3><ul>";
        o.carriers.forEach(function(t) {
            var r = o.transcript_records[t];
            h += "<li>" + esc(LAB.label(t)) + ": " + r.start + "–" + r.end + "; frame " + r.frame + "; " + esc(r.relationship) + "; CDS overlap " + r.cds_overlap_nt + " nt</li>";
        });
        h += "</ul>";
        h += "<h3>Protein annotations</h3>" + proteinDiagram(o);
        h += "<p>Physicochemical values are computational predictions, not measurements of cellular stability. Domain annotations are shown only for an exact sequence match to their cited database record.</p>";
        [ [ "Nucleotide sequence (including stop)", o.nt ], [ "Amino-acid sequence", o.aa ] ].forEach(function(pair) {
            h += '<div class="seqblock"><div class="seqhead"><span>' + esc(pair[0]) + '</span><button class="copy" data-seq="' + esc(pair[1]) + '">Copy</button></div><div class="seq">' + esc(pair[1]) + "</div></div>";
        });
        return h;
    }
    var SVGNS = "http://www.w3.org/2000/svg";
    function composeExport() {
        var GAP = 6, TITLEH = 48, totalH = TITLEH + IDEO_H + GAP + HDR_H + GAP + H;
        var out = document.createElementNS(SVGNS, "svg");
        out.setAttribute("width", W);
        out.setAttribute("height", totalH);
        out.setAttribute("viewBox", "0 0 " + W + " " + totalH);
        var style = document.createElementNS(SVGNS, "style");
        style.textContent = "text{font-family:" + FONT + "}";
        out.appendChild(style);
        var metadata = document.createElementNS(SVGNS, "metadata");
        metadata.textContent = JSON.stringify(TX_ORDER.map(function(t) {
            return {
                canonical_name: LAB.canonical(t),
                lab_alias: LAB.alias(t),
                accession: txAcc(t),
                source: txSrc(t),
                reference_release: txRel(t),
                assembly: META.assembly,
                exons: META.transcripts[t].exons,
                sequence: META.transcripts[t].sequence,
                cds: META.transcripts[t].cds,
                translation: META.transcripts[t].protein_sequence,
                evidence_status: META.transcripts[t].evidence_class,
                application_version: META.application_version,
                reference_build: META.reference_build
            };
        }));
        out.appendChild(metadata);
        var bg = document.createElementNS(SVGNS, "rect");
        bg.setAttribute("x", 0);
        bg.setAttribute("y", 0);
        bg.setAttribute("width", W);
        bg.setAttribute("height", totalH);
        bg.setAttribute("fill", "#ffffff");
        out.appendChild(bg);
        var title = document.createElementNS(SVGNS, "text");
        title.setAttribute("x", 14);
        title.setAttribute("y", 20);
        title.setAttribute("font-size", 12);
        title.textContent = "CURRENT · GENCODE " + META.gencode_release + " · " + META.assembly + " · RefSeq " + META.refseq_accession;
        out.appendChild(title);
        var version = document.createElementNS(SVGNS, "text");
        version.setAttribute("x", 14);
        version.setAttribute("y", 38);
        version.setAttribute("font-size", 10);
        version.textContent = "App " + META.application_version + " · " + META.reference_build;
        out.appendChild(version);
        [ [ "ideo", TITLEH ], [ "hdr", TITLEH + IDEO_H + GAP ], [ "stack", TITLEH + IDEO_H + GAP + HDR_H + GAP ] ].forEach(function(pair) {
            var src = $(pair[0]);
            if (!src) return;
            var g = document.createElementNS(SVGNS, "g");
            g.setAttribute("transform", "translate(0," + pair[1] + ")");
            Array.prototype.forEach.call(src.childNodes, function(n) {
                g.appendChild(n.cloneNode(true));
            });
            out.appendChild(g);
        });
        return {
            node: out,
            width: W,
            height: totalH
        };
    }
    function serializeExport() {
        var c = composeExport();
        return {
            xml: '<?xml version="1.0" encoding="UTF-8"?>\n' + (new XMLSerializer).serializeToString(c.node),
            width: c.width,
            height: c.height
        };
    }
    function triggerDownload(blob, name) {
        var url = URL.createObjectURL(blob), a = document.createElement("a");
        a.href = url;
        a.download = name;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        setTimeout(function() {
            URL.revokeObjectURL(url);
        }, 600);
    }
    function b64utf8(s) {
        return btoa(encodeURIComponent(s).replace(/%([0-9A-F]{2})/g, function(_, p) {
            return String.fromCharCode("0x" + p);
        }));
    }
    function downloadSVG() {
        try {
            triggerDownload(new Blob([ serializeExport().xml ], {
                type: "image/svg+xml;charset=utf-8"
            }), "MYCN_isoform_stack.svg");
        } catch (e) {
            console.error(e);
            toast("SVG export failed", true);
        }
    }
    function downloadPNG() {
        var ex, scale = 2.5;
        try {
            ex = serializeExport();
        } catch (e) {
            console.error(e);
            toast("PNG export failed", true);
            return;
        }
        var img = new Image;
        img.onload = function() {
            try {
                var canvas = document.createElement("canvas");
                canvas.width = Math.round(ex.width * scale);
                canvas.height = Math.round(ex.height * scale);
                var ctx = canvas.getContext && canvas.getContext("2d");
                if (!ctx) {
                    toast("PNG export failed (no canvas)", true);
                    return;
                }
                ctx.fillStyle = "#ffffff";
                ctx.fillRect(0, 0, canvas.width, canvas.height);
                ctx.setTransform(scale, 0, 0, scale, 0, 0);
                ctx.drawImage(img, 0, 0);
                if (!canvas.toBlob) {
                    toast("PNG export failed (no toBlob)", true);
                    return;
                }
                canvas.toBlob(function(blob) {
                    if (blob) triggerDownload(blob, "MYCN_isoform_stack.png"); else toast("PNG export failed", true);
                }, "image/png");
            } catch (e) {
                console.error(e);
                toast("PNG export failed", true);
            }
        };
        img.onerror = function() {
            toast("PNG export failed", true);
        };
        try {
            img.src = "data:image/svg+xml;base64," + b64utf8(ex.xml);
        } catch (e) {
            console.error(e);
            toast("PNG export failed", true);
        }
    }
    try {
        if (typeof globalThis !== "undefined") globalThis.__mstack = {
            setMETA: function(m, s, c) {
                META = m;
                LOCSEQ = s || LOCSEQ;
                CYTO = c || CYTO;
                if (META.palette) {
                    for (var k in META.palette) if (PAL[k] !== undefined) PAL[k] = META.palette[k];
                }
                TX_ORDER.forEach(function(t) {
                    state[t] = {
                        open: false,
                        sharedOpen: false
                    };
                });
                deriveIntrons();
                deriveGPos();
                buildLegendData();
            },
            acenSpan: function() {
                return acenSpan();
            },
            bandAt: function(p) {
                return bandAt(p);
            },
            STAIN: STAIN,
            verifyBaseMap: verifyBaseMap,
            baseAt: baseAt,
            specGroups: specGroups,
            exonSegments: exonSegments,
            orfIdsOfTx: orfIdsOfTx,
            orfsOfClass: orfsOfClass,
            carrierCount: carrierCount,
            txCoding: txCoding,
            txSpan: txSpan,
            renderCard: renderCard,
            badgesFor: badgesFor,
            codonGenomic: function(id, k) {
                var g = GPOS[id];
                return [ g[3 * k], g[3 * k + 1], g[3 * k + 2] ];
            },
            flyTo: function(a, b) {
                return flyTo(a, b, false);
            },
            view: function() {
                return viewG.slice();
            },
            ppb: pxPerBase,
            zoomToRange: function(a, b) {
                return zoomToRange(a, b, false);
            },
            brushEnd: function(px0, px1) {
                return onRulerBrushEnd({
                    selection: [ px0, px1 ]
                });
            },
            xOf: function(p) {
                return x(p);
            },
            invX: function(px) {
                return x.invert(px);
            },
            plotL: plotL,
            plotR: plotR,
            rightOpen: function() {
                return rightOpen;
            },
            openRight: openRight,
            noAnim: function(v) {
                NOANIM = !!v;
            },
            zoomFilter: zoomFilter,
            scrollStackBy: scrollStackBy,
            expandAll: function(full) {
                TX_ORDER.forEach(function(t) {
                    state[t].open = true;
                    if (full) state[t].sharedOpen = true;
                });
                rebuild();
                renderFeatureList();
            },
            pickOrf: function(id) {
                return pickOrf(id);
            },
            selectedOrf: function() {
                return selectedOrf;
            },
            serializeExport: function() {
                return serializeExport();
            },
            dmycnTx: dmycnTx,
            GPOS: function() {
                return GPOS;
            },
            RETAINED: function() {
                return RETAINED;
            },
            CANON_INTRONS: function() {
                return CANON_INTRONS;
            },
            LEGEND: function() {
                return LEGEND;
            },
            PAL: function() {
                return PAL;
            },
            TX_ORDER: TX_ORDER,
            BASECOL: BASECOL,
            GROUPS: GROUPS,
            NOVEL: NOVEL
        };
    } catch (e) {}
})();