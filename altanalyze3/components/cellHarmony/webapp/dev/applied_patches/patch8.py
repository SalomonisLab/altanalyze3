P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/templates/index.html"
src = open(P).read()

for panel in ("viz1", "viz2"):
    old = f'''                      <label class="field hidden inline-field compact-inline-field" id="{panel}-modality-field">
                        <span>Modality</span>
                        <select id="{panel}-modality" class="compact-select"></select>
                      </label>'''
    new = old + f'''
                      <!-- Color by. The UMAP cell-type view colours each cell by its
                           cellHarmony assignment. Any other categorical column of obs
                           works as well, so an atlas annotation the uploaded h5ad
                           already carries can be drawn instead. -->
                      <label class="field hidden inline-field compact-inline-field" id="{panel}-colorby-field">
                        <span>Color by</span>
                        <select id="{panel}-colorby" class="compact-select"></select>
                      </label>
                      <!-- Coordinates. Any 2-D embedding in obsm can replace the
                           cellHarmony projection. The reference atlas holds different
                           coordinates, so the app hides it while another embedding or
                           another colour column is selected. -->
                      <label class="field hidden inline-field compact-inline-field" id="{panel}-coords-field">
                        <span>Coordinates</span>
                        <select id="{panel}-coords" class="compact-select"></select>
                      </label>'''
    assert src.count(old) == 1, f"{panel}: {src.count(old)} matches"
    src = src.replace(old, new)

# --- chat examples come from the job, not from eight fixed lung sentences -----
old_examples = '''            <div class="chat-examples">
              <span>Try:</span>
              <button type="button" class="ghost-btn chat-example">Which genes are significant in COPD versus controls in AT2 cells?</button>
              <button type="button" class="ghost-btn chat-example">Is the AT2 COPD signature present in all donors or a subset?</button>
              <button type="button" class="ghost-btn chat-example">Which cell type is most affected in COPD?</button>
              <button type="button" class="ghost-btn chat-example">Which COPD genes change across many cell types rather than one?</button>
              <button type="button" class="ghost-btn chat-example">Do the clinical clusters differ molecularly in AT2?</button>
              <button type="button" class="ghost-btn chat-example">What distinguishes AT1 from AT2 cells?</button>
              <button type="button" class="ghost-btn chat-example">What are the best marker genes of AT2 cells?</button>
              <button type="button" class="ghost-btn chat-example">Show me the transcriptional targets of RUNX1 in alveolar macrophages</button>
            </div>'''
new_examples = '''            <!-- The examples name this job's own reference and cell states, so a
                 bone-marrow job never offers a lung question. app.js fills them from
                 /api/jobs/{id}/chat-examples when the results load. -->
            <div class="chat-examples" id="chat-examples">
              <span>Try:</span>
            </div>'''
assert src.count(old_examples) == 1
src = src.replace(old_examples, new_examples)

old_placeholder = '''              <textarea id="chat-question" rows="2"
                placeholder="e.g. What are the best marker genes of AT2 cells?"></textarea>'''
new_placeholder = '''              <textarea id="chat-question" rows="2"
                placeholder="Load a dataset to see example questions"></textarea>'''
assert src.count(old_placeholder) == 1
src = src.replace(old_placeholder, new_placeholder)

open(P, "w").write(src)
print("patch8 (index.html) applied")
