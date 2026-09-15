P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/app.py"
src = open(P).read()

anchor = '''    @app.api_route("/api/jobs/{job_id}/marker/heatmap.tsv", methods=["GET", "HEAD", "OPTIONS"])'''

routes = '''    @app.get("/api/jobs/{job_id}/chat-examples")
    def chat_examples(job_id: str):
        """Example questions for THIS job, named after its own cell states.

        The Chat tab shipped eight fixed lung sentences, so a bone-marrow job
        offered AT2 and COPD examples that its data cannot answer.
        """
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        return JSONResponse(_chat_examples(app, store.get_job(job_id)))

    @app.post("/api/jobs/{job_id}/chat")
    def chat(job_id: str, payload: ChatRequest):
        """Answer one question about this job, with numbers from this job.

        Two steps, deliberately separated. The assistant reads the sentence into
        one supported query and sees no data. This route then runs that query
        against the job's own files. So the model chooses the question and the
        data answers it; no number here comes from the model.
        """
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        question = str(payload.question or "").strip()
        if not question:
            raise HTTPException(status_code=400, detail="question is required")

        cache = _get_expression_cache(app, meta)
        reading = _chat_read_question(question, cache, meta)
        intent = str(reading.get("intent") or "")
        state = str(reading.get("cell_state") or "")
        state2 = str(reading.get("cell_state_2") or "")
        genes = [str(g) for g in (reading.get("genes") or []) if str(g).strip()]
        states = [s for s, _ in _chat_states_by_size(cache)]
        contrast = _chat_contrast(meta)
        result: Dict[str, Any] = {"question": question, "reading": reading, "intent": intent}

        if intent == "clarify":
            result["answer"] = (
                f"I need to know {reading.get('missing') or 'a little more'}. "
                f"This dataset has {len(states)} cell states and "
                f"{1 if contrast else 0} completed comparison(s).")
            result["choices"] = {"states": states,
                                 "contrasts": [contrast["label"]] if contrast else []}
            return JSONResponse(result)

        if intent == "unsupported":
            result["answer"] = (
                "I can answer four things about this dataset: the marker genes of a "
                "cell state, genes differing between two groups within a cell state, "
                "where named genes are expressed, and what separates two cell states.")
            return JSONResponse(result)

        if intent in _CHAT_NOT_YET:
            result["answer"] = (
                f"That is the {intent} protocol. It needs {_CHAT_NOT_YET[intent]}, which "
                "cellHarmony web does not compute yet, so I am not going to answer it "
                "with a different analysis.")
            result["status"] = "not_implemented"
            return JSONResponse(result)

        if intent in ("regulatory_driver", "communication_rewiring", "pathway_program"):
            label = {"regulatory_driver": "marker and differential networks",
                     "communication_rewiring": "cell communication",
                     "pathway_program": "GO Terms"}[intent]
            tab = "Explore" if intent == "communication_rewiring" else "Differential"
            result["answer"] = (
                f"That is the {intent} protocol, answered by this job's {label} results "
                f"for {state or 'the selected cell state'}. Open the {tab} tab for the "
                "figure; the chat does not inline it.")
            result["status"] = "use_existing_view"
            result["plot"] = {"kind": "network"}
            return JSONResponse(result)

        intent = _CHAT_PROTOCOL_ALIAS.get(intent, intent)

        if intent == "markers":
            if state and state not in states:
                result["answer"] = (f"{state} is not a cell state in this dataset. "
                                    f"It holds {len(states)}: {', '.join(states[:8])}"
                                    + (" ..." if len(states) > 8 else "."))
                result["status"] = "not_found"
                return JSONResponse(result)
            rows = _chat_markers_for_states(app, meta, cache, [state or states[0]], 25)
            sources = sorted({str(row["source"]) for row in rows})
            result["answer"] = (f"{len(rows)} top marker genes of {state or states[0]}, "
                                f"from the {' and '.join(sources) or 'marker'} analysis.")
            result["table"] = {
                "columns": ["gene", "cluster", "fold", "p", "source"],
                "rows": [[r["gene"], r["cluster"], round(float(r["fold"]), 4), r["p"], r["source"]]
                         for r in rows]}
            result["plot"] = {"kind": "dotplot", "genes": [r["gene"] for r in rows[:12]]}
            return JSONResponse(result)

        if intent == "compare":
            pair = [s for s in (state, state2) if s]
            if len(pair) < 2:
                result["answer"] = "Name two cell states to compare."
                result["status"] = "clarify"
                return JSONResponse(result)
            rows = _chat_markers_for_states(app, meta, cache, pair, 30)
            result["answer"] = f"Marker genes separating {pair[0]} and {pair[1]}."
            result["table"] = {
                "columns": ["gene", "cluster", "fold", "p", "source"],
                "rows": [[r["gene"], r["cluster"], round(float(r["fold"]), 4), r["p"], r["source"]]
                         for r in rows]}
            result["plot"] = {"kind": "dotplot", "genes": [r["gene"] for r in rows[:12]]}
            return JSONResponse(result)

        if intent == "expression":
            if not genes:
                result["answer"] = "Name at least one gene."
                result["status"] = "clarify"
                return JSONResponse(result)
            stats = _gene_state_stats(cache, genes)
            found = stats["genes"]
            if not found:
                result["answer"] = f"None of {', '.join(genes)} are in this dataset."
                result["status"] = "not_found"
                return JSONResponse(result)
            table_rows = []
            for row_index, gene in enumerate(found):
                means = np.asarray(stats["mean"][row_index], dtype=float)
                for column in np.argsort(-means)[:5]:
                    table_rows.append([gene, stats["states"][int(column)],
                                       round(float(means[int(column)]), 3),
                                       round(float(stats["frac"][row_index][int(column)]), 3)])
            missing = stats["missing"]
            result["answer"] = (
                f"Highest-expressing cell states for {', '.join(found)}, by mean expression "
                f"across the {len(stats['states'])} states of this dataset."
                + (f" Not in this dataset: {', '.join(missing)}." if missing else ""))
            result["table"] = {"columns": ["gene", "cell state", "mean", "fraction"],
                               "rows": table_rows}
            result["plot"] = {"kind": "dotplot", "genes": found}
            return JSONResponse(result)

        if intent == "differential":
            if not contrast:
                result["answer"] = (
                    "This job has not run a differential comparison yet. Open the "
                    "Differential tab, choose the two sample groups and run it; then this "
                    "question has an answer. I am not substituting a marker analysis for it.")
                result["status"] = "not_run"
                return JSONResponse(result)
            try:
                detail = _get_differential_detail_table(app, meta)
            except (FileNotFoundError, ValueError) as exc:
                result["answer"] = f"The differential result is unavailable ({exc})."
                result["status"] = "not_available"
                return JSONResponse(result)
            covered = sorted(set(detail["population"].astype(str)))
            subset = detail.loc[detail["population"].astype(str) == state] if state else detail
            if state and subset.empty:
                result["answer"] = (
                    f"{contrast['label']} is not computed for {state}. It covers "
                    f"{len(covered)} of the {len(states)} cell states in this dataset. "
                    "This is missing data, not an absence of change.")
                result["status"] = "not_covered"
                result["table"] = {"columns": ["cell state covered"], "rows": [[s] for s in covered]}
                return JSONResponse(result)
            frame = subset.copy()
            frame["fdr"] = pd.to_numeric(frame.get("fdr"), errors="coerce")
            frame["log2fc"] = pd.to_numeric(frame.get("log2fc"), errors="coerce")
            frame["pval"] = pd.to_numeric(frame.get("pval"), errors="coerce")
            frame = frame.dropna(subset=["log2fc"]).sort_values(["fdr", "pval"])
            top = frame.head(25)
            result["answer"] = (
                f"Top genes for {contrast['label']}" + (f" in {state}" if state else "")
                + f", from this job's cellHarmony-differential result "
                  f"({len(frame)} genes reported{' in ' + state if state else ''}).")
            result["table"] = {
                "columns": ["gene", "population", "log2fc", "fdr", "pval"],
                "rows": [[str(row.gene), str(row.population),
                          float(row.log2fc),
                          float(row.fdr) if _is_finite_number(row.fdr) else None,
                          float(row.pval) if _is_finite_number(row.pval) else None]
                         for row in top.itertuples()]}
            result["plot"] = {"kind": "volcano"}
            return JSONResponse(result)

        result["answer"] = "I did not understand that."
        return JSONResponse(result)

'''

assert src.count(anchor) == 1
src = src.replace(anchor, routes + anchor)
open(P, "w").write(src)
print("patch7 (chat routes) applied")
