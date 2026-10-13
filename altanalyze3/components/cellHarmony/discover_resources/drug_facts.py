"""What is known about a compound: FDA status, label warnings, adverse events, trials, papers.

Nathan, 2026-09-23: Interpret on Drug matching should say "what the relevance is for that
drug, in that cell type or disease context and if any of the drugs are FDA approved, if
there are adverse events or toxicity associated with any of the top ranking drugs or
which drugs the llm would suggest to proceed with given these variables".

IDENTITY FIRST. A compound is matched to the FDA by its UNII, the FDA's own substance
identifier, which PubChem lists among a compound's synonyms. The name alone failed: the
screen calls one compound fludroxycortide, the FDA files it as flurandrenolide (Cordran),
and a name search reported an approved drug as having no FDA application. Resolution is
the name first (PubChem's name lookup returned the right CID for all three test names),
then the screen's own PubChem CID, and the generic name only when neither gives a UNII.

  fda        openFDA drugsfda: every FDA application filed under that generic name, with
             its sponsor, brand names, marketing status and the earliest approval date
             openFDA lists, over up to 100 applications. Read over the first 20 only,
             cimetidine dated from 1994; over all 57 it dates from 1977 (Tagamet). A
             molecule whose original application openFDA does not list still dates late.
  label      openFDA drug label: the boxed warning, the indications, the mechanism.
  faers      openFDA adverse events (FAERS): the most reported MedDRA reactions and the
             number of reports. A report count is not an incidence: FAERS holds reports,
             not exposures. Product-use terms ("off label use", "drug ineffective") are
             returned in their own list, because they are not adverse reactions.
  trials     ClinicalTrials.gov API v2, the source LungChat's clinical_trials_search_tool
             uses, with the same query form: `query.intr` for the compound and
             `query.cond` for the disease. LungChat's ontology expansion reads an
             `ontologies.sqlite` its checkout does not ship, so the condition here is the
             disease name the atlas vocabulary gives.
  pubmed     the compound with the disease, and with lung, in title or abstract.

ChEMBL is not used. On 2026-09-23 its REST API answered HTTP 500 or timed out on every
probe, and a source that fails silently would read as "not approved".
"""
from __future__ import annotations

import json
import re
import urllib.error
import urllib.parse
import urllib.request
from functools import lru_cache
from typing import Any

OPENFDA = "https://api.fda.gov/drug/"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
CTGOV = "https://clinicaltrials.gov/api/v2/studies"
TIMEOUT = 10.0

# MedDRA terms that describe how a product was used, not what it did to the patient.
PRODUCT_USE = {
    "OFF LABEL USE", "DRUG INEFFECTIVE", "PRODUCT USE IN UNAPPROVED INDICATION",
    "PRODUCT USE ISSUE", "INTENTIONAL PRODUCT USE ISSUE", "NO ADVERSE EVENT",
    "INAPPROPRIATE SCHEDULE OF PRODUCT ADMINISTRATION", "PRODUCT DOSE OMISSION ISSUE",
    "PRODUCT DOSE OMISSION", "WRONG TECHNIQUE IN PRODUCT USAGE PROCESS",
    "INCORRECT DOSE ADMINISTERED", "DRUG DOSE OMISSION", "INTENTIONAL OVERDOSE",
    "PRODUCT QUALITY ISSUE", "EXPIRED PRODUCT ADMINISTERED", "THERAPY CESSATION",
    "THERAPY INTERRUPTED", "TREATMENT NONCOMPLIANCE", "DRUG INTERACTION",
    "THERAPEUTIC PRODUCT EFFECT INCOMPLETE", "THERAPEUTIC RESPONSE DECREASED",
    "OFF LABEL USE IN UNAPPROVED INDICATION", "CONDITION AGGRAVATED",
    "PRODUCT USE COMPLAINT", "DRUG EFFECTIVE FOR UNAPPROVED INDICATION",
}

# A screen identifier is not a name any registry files under.
_UNNAMED = re.compile(r"^(BRD|SA|NCGC|CID|CHEMBL|DB|LINCS|PRISM)[-_:]?\w*\d", re.I)


def nameable(name: str) -> bool:
    name = (name or "").strip()
    return bool(name) and not _UNNAMED.match(name) and bool(re.search(r"[A-Za-z]{4}", name))


def _json(url: str) -> tuple[dict | None, str]:
    """(payload, error). An openFDA 404 means no match and returns an empty payload."""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "LungMAP.net interpret"})
        with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
            return json.loads(r.read().decode()), ""
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return {}, ""
        return None, f"HTTP {exc.code}"
    except Exception as exc:                                           # noqa: BLE001
        return None, type(exc).__name__


def _clip(text: str, n: int) -> str:
    text = re.sub(r"\s+", " ", text or "").strip()
    return text if len(text) <= n else text[: n].rsplit(" ", 1)[0] + "…"


@lru_cache(maxsize=2048)
def identity(name: str, cid: str = "") -> dict[str, Any]:
    """The compound's PubChem CID and UNII, by name first and the screen's CID second."""
    tried = []
    for how, url in (("PubChem name", f"{PUBCHEM}name/{urllib.parse.quote(name)}/synonyms/JSON"),
                     ("PubChem CID", f"{PUBCHEM}cid/{cid}/synonyms/JSON" if str(cid).strip()
                      else "")):
        if not url:
            continue
        got, err = _json(url)
        tried.append(how + (f" ({err})" if err else ""))
        info = ((got or {}).get("InformationList") or {}).get("Information") or []
        for i in info[:1]:
            uniis = [x[5:] for x in i.get("Synonym") or [] if x.startswith("UNII-")]
            if uniis:
                return {"cid": i.get("CID"), "unii": uniis[:3], "via": how, "tried": tried}
    return {"cid": None, "unii": [], "via": "", "tried": tried}


def _fda_search(prefix: str, ident: dict, name: str, extra: tuple = ()) -> str:
    """One openFDA OR-search: the UNII, the generic name, and any extra name fields.

    THE UNII ALONE MISSES SALTS. PubChem gives the parent's UNII, and the FDA files
    toremifene as toremifene citrate (UNII 2498Y783QT, not 7NFE54O27T), so a UNII-only
    search reported three approved drugs as unapproved. The generic name matches the salt
    by word, and the active-ingredient name finds applications openFDA no longer
    annotates: cerivastatin's Baycol, withdrawn in 2001, is found only that way.
    """
    clauses = []
    if ident.get("unii"):
        clauses.append(f"{prefix}.unii:(" + "+".join(f'"{u}"' for u in ident["unii"]) + ")")
    clauses.append(f'{prefix}.generic_name:"{name}"')
    clauses.extend(f'{field}:"{name}"' for field in extra)
    return urllib.parse.quote("+".join(clauses), safe=":()+")


def _matched_by(results: list, ident: dict, name: str) -> str:
    """Which clause found the applications, so the reader can judge the match."""
    uniis = set(ident.get("unii") or [])
    low = name.lower()
    if any(uniis & set((a.get("openfda") or {}).get("unii") or []) for a in results):
        return f"UNII {', '.join(sorted(uniis))} via {ident['via']}"
    if any(low in " ".join((a.get("openfda") or {}).get("generic_name") or []).lower()
           for a in results):
        return "generic name (salt forms included)"
    if results:
        return "active-ingredient name (salt forms and unannotated applications included)"
    return (f"UNII {', '.join(sorted(uniis))} via {ident['via']}, generic name and "
            "active-ingredient name" if uniis else "generic and active-ingredient name; "
            "PubChem gave no UNII")


@lru_cache(maxsize=2048)
def fda(name: str, cid: str = "") -> dict[str, Any]:
    """FDA applications, label and FAERS reactions for this compound, matched by UNII."""
    ident = identity(name, cid)
    qa = _fda_search("openfda", ident, name, ("products.active_ingredients.name",))
    q = _fda_search("openfda", ident, name)
    qe = _fda_search("patient.drug.openfda", ident, name, ("patient.drug.medicinalproduct",))
    out: dict[str, Any] = {"name": name, "errors": [],
                           "cid": ident.get("cid"), "unii": ident.get("unii")}

    # 100 applications, not 20: the earliest approval is a minimum over them, and
    # etodolac has 42, so a first page of 20 could have missed the 1991 original.
    apps, err = _json(f"{OPENFDA}drugsfda.json?search={qa}&limit=100")
    out["matched_by"] = _matched_by((apps or {}).get("results") or [], ident, name)
    if err:
        out["errors"].append(f"drugsfda {err}")
    applications = []
    for a in (apps or {}).get("results") or []:
        dates = [s.get("submission_status_date") for s in a.get("submissions") or []
                 if s.get("submission_type") == "ORIG" and s.get("submission_status") == "AP"
                 and s.get("submission_status_date")]
        products = a.get("products") or []
        applications.append({
            "application": a.get("application_number", ""),
            "sponsor": a.get("sponsor_name", ""),
            "brands": sorted({p.get("brand_name", "") for p in products if p.get("brand_name")}),
            "status": sorted({p.get("marketing_status", "") for p in products
                              if p.get("marketing_status")}),
            "approved": min(dates) if dates else "",
        })
    out["applications"] = applications
    out["n_applications_read"] = len(applications)
    out["n_applications"] = int(((apps or {}).get("meta") or {}).get("results", {})
                                .get("total", len(applications)) or 0)
    dated = sorted(a["approved"] for a in applications if a["approved"])
    out["earliest_listed_approval"] = (f"{dated[0][:4]}-{dated[0][4:6]}-{dated[0][6:8]}"
                                       if dated else "")
    out["fda_names"] = sorted({i["name"].lower() for a in (apps or {}).get("results") or []
                               for p in a.get("products") or []
                               for i in p.get("active_ingredients") or [] if i.get("name")}
                              | {g.lower() for a in (apps or {}).get("results") or []
                                 for g in (a.get("openfda") or {}).get("generic_name") or []})
    # NDA and BLA are new-drug applications; ANDA is a generic of an approved drug. Any of
    # the three with an approval date is an FDA approval of this active ingredient.
    out["fda_approved"] = bool(applications) if apps is not None else None
    out["marketed"] = sorted({s for a in applications for s in a["status"]})
    out["url"] = ("https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event="
                  "BasicSearch.process&searchterm=" + urllib.parse.quote(name))

    # EVERY LABEL, NOT THE FIRST. One label per manufacturer: etodolac has 73, and a label
    # taken at random came back without the boxed-warning field on one call and with it on
    # the next, which moved the drug between groups. The boxed warning is counted over the
    # labels read and quoted from the newest label that carries one.
    lab, err = _json(f"{OPENFDA}label.json?search={q}&limit=100")
    if err:
        out["errors"].append(f"label {err}")
    labels = (lab or {}).get("results") or []
    out["label_found"] = bool(labels)
    out["n_labels"] = int(((lab or {}).get("meta") or {}).get("results", {})
                          .get("total", len(labels)) or 0)
    out["n_labels_read"] = len(labels)
    boxed = sorted((x for x in labels if x.get("boxed_warning")),
                   key=lambda x: x.get("effective_time") or "", reverse=True)
    out["n_labels_boxed"] = len(boxed)
    newest = sorted(labels, key=lambda x: x.get("effective_time") or "", reverse=True)
    r = boxed[0] if boxed else (newest[0] if newest else {})
    out["boxed_warning"] = _clip(" ".join(r.get("boxed_warning") or []), 420)
    # A label section opens with its own number and heading ("12.1 Mechanism of Action").
    head = r"^\s*[\d.]*\s*(mechanism of action|indications and usage|indications & usage)\s*:?\s*"
    out["indications"] = _clip(re.sub(head, "", " ".join(r.get("indications_and_usage") or []),
                                      flags=re.I), 300)
    out["mechanism"] = _clip(re.sub(head, "", " ".join(r.get("mechanism_of_action") or []),
                                    flags=re.I), 260)
    out["label_url"] = (f"https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid={r['set_id']}"
                        if r.get("set_id") else "")

    ev, err = _json(f"{OPENFDA}event.json?search={qe}"
                    "&count=patient.reaction.reactionmeddrapt.exact&limit=25")
    tot, err2 = _json(f"{OPENFDA}event.json?search={qe}&limit=1")
    for e in (err, err2):
        if e:
            out["errors"].append(f"FAERS {e}")
    terms = (ev or {}).get("results") or []
    out["faers_reports"] = int(((tot or {}).get("meta") or {}).get("results", {})
                               .get("total", 0) or 0)
    out["faers_reactions"] = [{"term": t["term"].capitalize(), "count": t["count"]}
                              for t in terms if t["term"] not in PRODUCT_USE][:6]
    out["faers_product_use"] = [{"term": t["term"].capitalize(), "count": t["count"]}
                                for t in terms if t["term"] in PRODUCT_USE][:4]
    out["faers_url"] = "https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html"
    return out


@lru_cache(maxsize=2048)
def trials(name: str, condition: str = "") -> dict[str, Any]:
    """ClinicalTrials.gov studies of this intervention, and of it in this condition.

    `name` may be an OR of the screen's name and the FDA's ("fludroxycortide OR
    flurandrenolide"), which ClinicalTrials.gov's Essie syntax accepts.
    """
    out: dict[str, Any] = {"name": name, "condition": condition, "errors": []}
    base = {"fields": "NCTId,BriefTitle,OverallStatus,Phase,Condition",
            "countTotal": "true", "format": "json"}
    every, err = _json(CTGOV + "?" + urllib.parse.urlencode(
        dict(base, **{"query.intr": name, "pageSize": 100})))
    if err:
        out["errors"].append(f"ClinicalTrials.gov {err}")
    rank = {"EARLY_PHASE1": 0.5, "PHASE1": 1, "PHASE2": 2, "PHASE3": 3, "PHASE4": 4}
    phases = [p for s in (every or {}).get("studies") or []
              for p in (s["protocolSection"].get("designModule") or {}).get("phases") or []
              if p in rank]
    out["n_trials"] = int((every or {}).get("totalCount") or 0)
    out["n_read"] = len((every or {}).get("studies") or [])
    out["max_phase"] = (max(phases, key=rank.get).replace("_", " ").title()
                        .replace("Phase", "Phase ") if phases else "")
    out["url"] = ("https://clinicaltrials.gov/search?" +
                  urllib.parse.urlencode({"intr": name}))
    if condition:
        hit, err = _json(CTGOV + "?" + urllib.parse.urlencode(
            dict(base, **{"query.intr": name, "query.cond": condition, "pageSize": 3})))
        if err:
            out["errors"].append(f"ClinicalTrials.gov {err}")
        out["n_in_condition"] = int((hit or {}).get("totalCount") or 0)
        out["in_condition"] = [
            {"nct": s["protocolSection"]["identificationModule"]["nctId"],
             "title": s["protocolSection"]["identificationModule"].get("briefTitle", ""),
             "phase": ", ".join((s["protocolSection"].get("designModule") or {})
                                .get("phases") or []),
             "status": (s["protocolSection"].get("statusModule") or {})
             .get("overallStatus", "")}
            for s in (hit or {}).get("studies") or []]
        out["condition_url"] = ("https://clinicaltrials.gov/search?" +
                                urllib.parse.urlencode({"intr": name, "cond": condition}))
    return out


def pubmed(name: str, disease: str = "") -> dict[str, Any]:
    """The cached lookup; a failure is returned, not cached."""
    try:
        return _pubmed(name, disease)
    except Exception as exc:                                           # noqa: BLE001
        code = getattr(exc, "code", "")
        return {"themes": {}, "error": f"PubMed could not be reached "
                f"({type(exc).__name__}{' ' + str(code) if code else ''})."}


@lru_cache(maxsize=2048)
def _pubmed(name: str, disease: str = "") -> dict[str, Any]:
    """PubMed counts for the compound with the disease and with lung."""
    from .gene_literature import EUTILS, _get, phrase_terms, pubmed_url

    names = [re.sub(r"[^A-Za-z0-9 -]", " ", n).strip() for n in name.split(" OR ")]
    drug = " OR ".join(f'"{n}"[tiab]' for n in names if n)
    themes = []
    if disease:
        themes.append(("disease", f"({drug}) AND ({phrase_terms(disease)})"))
    themes.append(("lung", f"({drug}) AND (lung[tiab] OR pulmonary[tiab])"))
    out: dict[str, Any] = {"themes": {}}
    for theme, query in themes:
        got = _get(f"{EUTILS}esearch.fcgi?db=pubmed&retmode=json&retmax=0&term="
                   + urllib.parse.quote(query))["esearchresult"]
        out["themes"][theme] = {"count": int(got.get("count") or 0), "query": query,
                                "url": pubmed_url(query)}
    return out


def facts(name: str, disease: str = "", cid: str = "") -> dict[str, Any]:
    if not nameable(name):
        return {"name": name, "unnamed": True,
                "note": f"{name} is a screen identifier, not a drug name, so no registry "
                        "files it."}
    fd = fda(name, str(cid or ""))
    # The FDA's own generic names join the screen's for the trial and paper searches.
    # A salt ("toremifene citrate") adds nothing to a word search for "toremifene"; a
    # different name ("flurandrenolide" for fludroxycortide) does.
    names = [name] + [n for n in fd.get("fda_names") or []
                      if name.lower() not in n and len(n.split()) <= 2][:2]
    both = " OR ".join(names)
    return {"name": name, "fda": fd, "searched_as": names,
            "trials": trials(both, disease), "pubmed": pubmed(both, disease)}


# Drug names end in these stems (USAN/INN). A recommendation that names such a word not in
# the facts is naming a drug from memory, and is withheld.
_DRUG_STEM = re.compile(
    r"\b[a-z]+(?:mab|nib|lib|vir|mycin|micin|cillin|olol|pril|sartan|azole|afil|"
    r"parin|dronate|lukast|oxacin|cycline|tinib|zumab|ximab|umab|platin|taxel|"
    r"rubicin|profen|coxib|semide|thiazide|dipine|triptan|setron|prazole|gliptin|"
    r"glitazone|fibrate|olone|asone|zepam|barbital)\b", re.I)
# Stems that also end ordinary words ("-terol" in cholesterol, "-tide" in peptide, "-stat"
# in thermostat) are left out, so the guard does not withhold a correct sentence.


_RECOMMENDED: dict[str, dict] = {}


def recommend(setting: str, question: str, rows: list[dict]) -> dict[str, Any]:
    """The local model's pick among the compounds on screen, from the facts only.

    Kept per setting, question and compound list, because generation takes about 45 s.
    """
    key = json.dumps([setting, question, [(r["compound"], r["best_rank"]) for r in rows]])
    if key not in _RECOMMENDED:
        if len(_RECOMMENDED) > 256:
            _RECOMMENDED.clear()
        _RECOMMENDED[key] = _recommend(setting, question, rows)
    return _RECOMMENDED[key]


def _fact_line(r: dict) -> str:
    """One compound's facts, self-contained, negatives stated, for the model and the guard."""
    f = r.get("facts") or {}
    fd, tr, pm = f.get("fda") or {}, f.get("trials") or {}, f.get("pubmed") or {}
    c = r["compound"]
    bits = [f"{c}: screen rank {r['best_rank']}, Spearman {r['rho']:+.2f} ({r['context']})"]
    if fd.get("fda_approved"):
        status = fd.get("marketed") or []
        bits.append("FDA approved" + ("; every FDA application discontinued"
                                      if status == ["Discontinued"] else "; marketed"))
    elif fd.get("fda_approved") is False:
        bits.append("no FDA approval found")
    bits.append(f"boxed warning: {_clip(fd['boxed_warning'], 80)}" if fd.get("boxed_warning")
                else "no boxed warning")
    if fd.get("faers_reactions"):
        bits.append("most reported reactions: " + ", ".join(
            t["term"].lower() for t in fd["faers_reactions"][:3]))
    bits.append(f"{tr['n_trials']} registered trials, highest {tr.get('max_phase') or 'phase unstated'}, "
                f"{tr.get('n_in_condition', 0)} in this disease" if tr.get("n_trials")
                else "no registered trials")
    d = (pm.get("themes") or {}).get("disease")
    if d is not None:
        bits.append(f"{d['count']} papers with this disease")
    return "; ".join(bits) + "."


_FUNCTION_WORDS = set("""with this that these those from have been also which their into
than more most such only each other both while when where there because being would could
should might about among over under after before within without""".split())


def _cross_attributed(text: str, rows: list[dict], lines: dict[str, str],
                      given: str) -> list[str]:
    """Sentences that tie a compound to a fact found only in another compound's line.

    The local model moved toremifene's QT boxed warning onto cerivastatin. A sentence
    naming one compound may use only words from that compound's own line, the setting,
    the question, or none of the lines at all; a word that occurs only in OTHER lines is a
    misattribution.
    """
    bad = []
    names = [r["compound"].lower() for r in rows if r["compound"].lower() in lines]
    for sentence in re.split(r"(?<=[.;])\s+", text or ""):
        low = sentence.lower()
        hit = [n for n in names if n in low]
        if len(hit) != 1:
            continue
        own = lines[hit[0]]
        others = " ".join(v for k, v in lines.items() if k != hit[0])
        for w in set(re.findall(r"[a-z][a-z-]{3,}", low)) - _FUNCTION_WORDS:
            if w in own or w in given or w not in others:
                continue
            bad.append(f"{hit[0]}: \"{w}\"")
    return sorted(set(bad))


def _recommend(setting: str, question: str, rows: list[dict]) -> dict[str, Any]:
    from .assistant import model

    if not model.available():
        return {"text": "", "note": "The local model is not loaded."}
    usable = [r for r in rows if not (r.get("facts") or {}).get("unnamed")]
    lines = {r["compound"].lower(): _fact_line(r).lower() for r in usable}
    context = "\n".join([setting] + [_fact_line(r) for r in usable])
    text = model.answer(question, context, max_tokens=130)
    # A trailing list number ("... 2.") is a cut-off list item, not a sentence.
    text = re.sub(r"\s*\d+\.\s*$", "", text or "").strip()
    # The question is part of what the model was given: "AT2" is in it, not in the facts.
    given = (setting + "\n" + question).lower()
    allowed = (context + "\n" + question)
    low = allowed.lower()
    stray = sorted({w for w in _DRUG_STEM.findall(text) if w.lower() not in low})
    upper = sorted({w for w in re.findall(r"\b[A-Z][A-Z0-9-]{2,}\b", text)
                    if w not in allowed})
    named = [r["compound"] for r in usable if r["compound"].lower() in text.lower()]
    crossed = _cross_attributed(text, usable, lines, given)
    if not text or stray or upper or not named or crossed:
        if crossed:
            why = ("it attributed a fact to the wrong compound (" + "; ".join(crossed[:3])
                   + " occurs only in another compound's facts)")
        elif stray or upper:
            why = (f"it named {', '.join((stray + upper)[:3])}, which the facts above do "
                   "not contain")
        else:
            why = "it named none of the compounds shown"
        return {"text": "", "withheld": text,
                "note": f"The model's suggestion was withheld because {why}."}
    return {"text": text, "named": named,
            "note": ("Written by the local model from the screen ranks and the FDA, "
                     "FAERS, ClinicalTrials.gov and PubMed facts above, and checked so that "
                     "no sentence ties a compound to another compound's facts. A "
                     "hypothesis for discussion, not a recommendation to treat.")}


def groups(rows: list[dict]) -> list[dict]:
    """The compounds grouped by the variables the reader weighs, in screen-rank order.

    A grouping, not a score: each group is a fact read from openFDA, and no weight is put
    on one variable against another.
    """
    order = [("approved_clean", "FDA approved and marketed, no boxed warning"),
             ("approved_boxed", "FDA approved and marketed, with a boxed warning"),
             ("discontinued", "FDA approved, every application now discontinued"),
             ("not_approved", "No FDA approval found"),
             ("unknown", "openFDA did not answer")]
    found: dict[str, list[str]] = {k: [] for k, _ in order}
    for r in rows:
        fd = (r.get("facts") or {}).get("fda") or {}
        if fd.get("fda_approved") is None:
            key = "unknown"
        elif not fd["fda_approved"]:
            key = "not_approved"
        elif (fd.get("marketed") or []) == ["Discontinued"]:
            key = "discontinued"
        elif fd.get("boxed_warning"):
            key = "approved_boxed"
        else:
            key = "approved_clean"
        found[key].append(r["compound"])
    return [{"key": k, "label": label, "compounds": found[k]} for k, label in order
            if found[k]]


_ONE: dict[str, dict] = {}


def assess_one(disease: str, cells: str, row: dict) -> dict[str, Any]:
    """The local model's reading of ONE compound, from that compound's facts alone.

    Asked to choose among eight compounds at once, the 1.5B model either copied the facts
    back or moved toremifene's QT boxed warning onto cerivastatin. Given one compound's
    line, it has nothing to confuse it with.
    """
    from .assistant import model

    line = _fact_line(row)
    key = json.dumps([disease, cells, line])
    if key in _ONE:
        return _ONE[key]
    c = row["compound"]
    if not model.available():
        return {"text": "", "note": "The local model is not loaded."}
    question = (f"In one or two sentences: is {c} worth testing against {disease} in "
                f"{cells} cells, and what is the main caution in its facts?")
    if len(question) > 200:
        question = (f"In one or two sentences: is {c} worth testing here, and what is the "
                    "main caution in its facts?")
    context = (f"Disease: {disease}. Cells: {cells}. A negative Spearman means the "
               f"compound's signature opposes the disease program.\n{line}")
    text = model.answer(question, context, max_tokens=70)
    text = re.sub(r"\s*\d+\.\s*$", "", text or "").strip()
    allowed = context + "\n" + question
    stray = sorted({w for w in _DRUG_STEM.findall(text) if w.lower() not in allowed.lower()})
    upper = sorted({w for w in re.findall(r"\b[A-Z][A-Z0-9-]{2,}\b", text)
                    if w not in allowed})
    echo = text and (text.lower().startswith(line.lower()[:40]) or line.lower()[:60] in
                     text.lower())
    if not text or stray or upper or echo or c.lower() not in text.lower():
        why = ("it copied the facts back" if echo else
               f"it named {', '.join((stray + upper)[:3])}, which the facts do not "
               "contain" if (stray or upper) else
               "it did not name the compound" if text else "it wrote nothing")
        out = {"text": "", "note": f"Withheld: {why}."}
    else:
        out = {"text": text, "note": ""}
    if len(_ONE) > 1024:
        _ONE.clear()
    _ONE[key] = out
    return out
