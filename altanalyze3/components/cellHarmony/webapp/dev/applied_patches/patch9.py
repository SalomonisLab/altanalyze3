P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/static/app.js"
src = open(P).read()

def sub(old, new):
    global src
    assert src.count(old) == 1, f"{src.count(old)} matches for:\n{old[:200]}"
    src = src.replace(old, new)

# --- 1. show the two new controls on the UMAP cell-type view ------------------
sub('''    if (combMinField) combMinField.classList.toggle("hidden", mode !== "combplot");
    if (wantsGeneSet) geneField.classList.add("hidden");''',
    '''    if (combMinField) combMinField.classList.toggle("hidden", mode !== "combplot");
    if (wantsGeneSet) geneField.classList.add("hidden");

    // The UMAP cell-type view colours cells by the cellHarmony assignment and
    // draws them on the cellHarmony projection. Both may be swapped for any
    // categorical obs column and any 2-D embedding the h5ad carries. The other
    // plot types have no such choice, so the two lists stay hidden there.
    const wantsUmapOptions = mode === "cluster";
    const colorByField = document.getElementById(panelElementId(panelKey, "colorby-field"));
    const coordsField = document.getElementById(panelElementId(panelKey, "coords-field"));
    if (colorByField) colorByField.classList.toggle("hidden", !wantsUmapOptions);
    if (coordsField) coordsField.classList.toggle("hidden", !wantsUmapOptions);
    if (wantsUmapOptions) refreshUmapOptions(panelKey);''')

# --- 2. the option lists ------------------------------------------------------
sub('''function appendGeneSetGroupParams(params, panelKey) {''',
    '''/* The Color by and Coordinates lists for the UMAP cell-type view.
 *
 * Both come from /plot-variables, which offers the categorical obs columns with
 * a workable number of levels and every obsm entry with two or more columns. The
 * empty value is the default in each list: the cellHarmony cell-state assignment
 * and the cellHarmony projection, which is what the panel drew before.
 */
async function refreshUmapOptions(panelKey) {
  const jobId = getResultsJobId();
  const colorBy = document.getElementById(panelElementId(panelKey, "colorby"));
  const coords = document.getElementById(panelElementId(panelKey, "coords"));
  if (!jobId || !colorBy || !coords) return;
  if (colorBy.dataset.jobId === jobId && coords.dataset.jobId === jobId) return;
  let info;
  try {
    info = await loadPlotVariables(jobId);
  } catch (err) {
    // Say it rather than leaving two empty lists that look broken.
    console.warn("UMAP options unavailable:", err);
    colorBy.innerHTML = '<option value="">unavailable</option>';
    coords.innerHTML = '<option value="">unavailable</option>';
    return;
  }

  const previousColor = colorBy.value;
  colorBy.innerHTML = "";
  const clusterOption = document.createElement("option");
  clusterOption.value = "";
  clusterOption.textContent = `${info.cluster_key} (cellHarmony)`;
  colorBy.appendChild(clusterOption);
  (info.color_variables || info.variables || []).forEach((variable) => {
    if (String(variable.field) === String(info.cluster_key)) return;
    const option = document.createElement("option");
    option.value = variable.field;
    option.textContent = `${variable.field} (${variable.n})`;
    colorBy.appendChild(option);
  });
  colorBy.value = Array.from(colorBy.options).some((o) => o.value === previousColor) ? previousColor : "";
  colorBy.dataset.jobId = jobId;

  const previousCoords = coords.value;
  coords.innerHTML = "";
  const coordOptions = (info.coords || []).length ? info.coords : [{ key: "", label: "cellHarmony UMAP" }];
  coordOptions.forEach((entry) => {
    const option = document.createElement("option");
    option.value = entry.key || "";
    option.textContent = entry.label || entry.key || "cellHarmony UMAP";
    coords.appendChild(option);
  });
  coords.value = Array.from(coords.options).some((o) => o.value === previousCoords) ? previousCoords : "";
  coords.dataset.jobId = jobId;

  if (colorBy.dataset.wired !== "1") {
    colorBy.dataset.wired = "1";
    colorBy.addEventListener("change", () => loadVisualizationPanel(panelKey));
  }
  if (coords.dataset.wired !== "1") {
    coords.dataset.wired = "1";
    coords.addEventListener("change", () => loadVisualizationPanel(panelKey));
  }
}

function appendGeneSetGroupParams(params, panelKey) {''')

# --- 3. the request carries the two choices ----------------------------------
sub('''    if (isUmapMode(mode)) {
      const params = getDisplayFilterParams(panelKey);
      params.set("modality", modality);''',
    '''    if (isUmapMode(mode)) {
      const params = getDisplayFilterParams(panelKey);
      params.set("modality", modality);
      // Only the cell-type view offers these. "UMAP broad" draws reference
      // against query, and the reference is hidden whenever either choice is
      // off-default, which would leave that view with nothing to contrast.
      if (mode === "cluster") {
        const colorBy = getPanelSelectValue(panelKey, "colorby");
        const coordsKey = getPanelSelectValue(panelKey, "coords");
        if (colorBy) params.set("color_by", colorBy);
        if (coordsKey) params.set("coords", coordsKey);
      }''')

# --- 4. say what is drawn, and why the reference is not there -----------------
sub('''  const populations = buildStableUmapPopulationOrder(umapData);
  const colorMap = buildReferencePreviewColorMap(populations);''',
    '''  setUmapPanelSummary(panelKey, umapData);
  const populations = buildStableUmapPopulationOrder(umapData);
  const colorMap = buildReferencePreviewColorMap(populations);''')

sub('''function renderPanelUmap(panelKey, umapData, mode, dotScale) {''',
    '''/* The caption under a UMAP panel. It names the colour column and the embedding
 * whenever they are not the cellHarmony defaults, and says why the reference
 * atlas is missing, which otherwise reads as a drawing fault. */
function setUmapPanelSummary(panelKey, umapData) {
  const parts = [];
  const filters = getDisplayFilterSummary(panelKey);
  if (filters) parts.push(`Display only: ${filters}`);
  if (umapData && umapData.reference_hidden) {
    parts.push(`Colored by ${umapData.color_label}`);
    parts.push(`Coordinates: ${umapData.coords_label}`);
    parts.push("Reference atlas hidden: it carries neither this annotation nor these coordinates");
  }
  if (parts.length) setPanelSummary(panelKey, parts.join(" | "));
}

function renderPanelUmap(panelKey, umapData, mode, dotScale) {''')

# --- 5. chat examples for the job that is loaded ------------------------------
sub('''    exploreResultsReadyJobId = normalizedJobId;
    exploreResultsPendingJobId = null;''',
    '''    exploreResultsReadyJobId = normalizedJobId;
    exploreResultsPendingJobId = null;
    void loadChatExamples(normalizedJobId);''')

sub('''async function askChat() {''',
    '''/* The Chat examples, built by the server from this job's reference and its own
 * cell states. The tab shipped eight fixed lung sentences, so a bone-marrow job
 * offered AT2 and COPD questions that its data cannot answer. */
async function loadChatExamples(jobId) {
  const host = document.getElementById("chat-examples");
  const box = document.getElementById("chat-question");
  if (!host || !jobId) return;
  try {
    const response = await fetch(apiPath(`/jobs/${jobId}/chat-examples`));
    const data = await response.json();
    if (!response.ok) throw new Error(data.detail || "examples unavailable");
    host.innerHTML = "";
    const label = document.createElement("span");
    label.textContent = data.reference ? `Try (${data.reference}):` : "Try:";
    host.appendChild(label);
    (data.examples || []).forEach((question) => {
      const button = document.createElement("button");
      button.type = "button";
      button.className = "ghost-btn chat-example";
      button.textContent = question;
      button.addEventListener("click", () => {
        box.value = question;
        askChat();
      });
      host.appendChild(button);
    });
    if (box && data.placeholder) box.placeholder = data.placeholder;
  } catch (err) {
    console.warn("chat examples unavailable:", err);
  }
}

async function askChat() {''')

open(P, "w").write(src)
print("patch9 (app.js) applied")
