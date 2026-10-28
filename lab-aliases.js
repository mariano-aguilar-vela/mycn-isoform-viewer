/* Owner-supplied laboratory terminology, 2026-10-28. Database identifiers remain authoritative. */
(function(root) {
  "use strict";
  var aliases = Object.freeze({"MYCN-201":"TV2", "MYCN-202":"TV3", "MYCN-203":"NC", "MYCN-204":"TV6", "MYCN-205":"TV1", "MYCN-206":"TV5", "MYCN-207":"TV2+", "MYCN-XM":"TV4"});
  function canonical(id) { return id.indexOf("MYCN-") === 0 ? id : "MYCN-" + id; }
  function alias(id) { return aliases[canonical(id)] || ""; }
  function label(id) { return canonical(id) + (alias(id) ? " (" + alias(id) + ")" : ""); }
  function accessible(id) { return canonical(id) + (alias(id) ? ", laboratory alias " + alias(id) : ""); }
  function matches(id, accession, query, extra) {
    var q = String(query || "").trim().toLowerCase();
    if (!q) return true;
    // Laboratory tokens are whole labels. In particular TV2 must never select TV2+.
    if (/^(?:tv\d+\+?|nc)$/.test(q)) return alias(id).toLowerCase() === q;
    return [canonical(id), accession, extra || ""].some(function(s) { return String(s).toLowerCase().indexOf(q) >= 0; });
  }
  var api = Object.freeze({aliases:aliases, canonical:canonical, alias:alias, label:label, accessible:accessible, matches:matches});
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.MYCNLab = api;
})(typeof window !== "undefined" ? window : globalThis);
