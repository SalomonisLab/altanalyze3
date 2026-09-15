P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/static/app.js"
src = open(P).read()

def sub(old, new):
    global src
    assert src.count(old) == 1, f"{src.count(old)} matches for:\n{old[:200]}"
    src = src.replace(old, new)

# --- 1. the coordinate list gains the obs-columns entry, and the axis lists ---
sub('''  const previousCoords = coords.value;
  coords.innerHTML = "";
  const coordOptions = (info.coords || []).length ? info.coords : [{ key: "", label: "cellHarmony UMAP" }];
  coordOptions.forEach((entry) => {
    const option = document.createElement("option");
    option.value = entry.key || "";
    option.textContent = entry.label || entry.key || "cellHarmony UMAP";
    coords.appendChild(option);
  });
  coords.value = Array.from(coords.options).some((o) => o.value === previousCoords) ? previousCoords : "";
  coords.dataset.jobId = jobId;''',
    '''  const previousCoords = coords.value;
  coords.innerHTML = "";
  const coordOptions = (info.coords || []).length ? info.coords : [{ key: "", label: "cellHarmony UMAP" }];
  const numeric = info.numeric_variables || [];
  // Two numeric annotations are the minimum for a pair of axes.
  const coordEntries = numeric.length > 1
    ? coordOptions.concat([{ key: OBS_AXES_KEY, label: "obs columns (pick X and Y)" }])
    : coordOptions;
  coordEntries.forEach((entry) => {
    const option = document.createElement("option");
    option.value = entry.key || "";
    option.textContent = entry.label || entry.key || "cellHarmony UMAP";
    coords.appendChild(option);
  });
  coords.value = Array.from(coords.options).some((o) => o.value === previousCoords) ? previousCoords : "";
  coords.dataset.jobId = jobId;

  // The axis lists. A field shows its range and how many cells carry a value,
  // because a metadata column recorded for part of the dataset draws a panel
  // with fewer cells than the UMAP has.
  const axisSelects = [
    [document.getElementById(panelElementId(panelKey, "xfield")), 0],
    [document.getElementById(panelElementId(panelKey, "yfield")), 1],
  ];
  axisSelects.forEach(([select, defaultIndex]) => {
    if (!select) return;
    const previous = select.value;
    select.innerHTML = "";
    numeric.forEach((entry) => {
      const option = document.createElement("option");
      option.value = entry.field;
      const missing = Number(entry.n_missing) > 0 ? `, ${entry.n_missing} without a value` : "";
      option.textContent = `${entry.field} (${formatAxisNumber(entry.min)} to ${formatAxisNumber(entry.max)}${missing})`;
      select.appendChild(option);
    });
    if (Array.from(select.options).some((o) => o.value === previous)) {
      select.value = previous;
    } else if (numeric[defaultIndex]) {
      select.value = numeric[defaultIndex].field;
    }
    select.dataset.jobId = jobId;
    if (select.dataset.wired !== "1") {
      select.dataset.wired = "1";
      select.addEventListener("change", () => loadVisualizationPanel(panelKey));
    }
  });''')

sub('''  if (coords.dataset.wired !== "1") {
    coords.dataset.wired = "1";
    coords.addEventListener("change", () => loadVisualizationPanel(panelKey));
  }
}''',
    '''  if (coords.dataset.wired !== "1") {
    coords.dataset.wired = "1";
    coords.addEventListener("change", () => {
      syncUmapAxisFields(panelKey);
      loadVisualizationPanel(panelKey);
    });
  }
  syncUmapAxisFields(panelKey);
}

/* The X and Y lists belong to the "obs columns" coordinate choice alone, so they
 * stay hidden while the panel draws an embedding. */
function syncUmapAxisFields(panelKey) {
  const mode = getPanelSelectValue(panelKey, "mode");
  const usingObsAxes = mode === "cluster" && getPanelSelectValue(panelKey, "coords") === OBS_AXES_KEY;
  ["xfield-field", "yfield-field"].forEach((suffix) => {
    const field = document.getElementById(panelElementId(panelKey, suffix));
    if (field) field.classList.toggle("hidden", !usingObsAxes);
  });
}

/* Axis ranges read as numbers a person can compare, not as 17 decimal places. */
function formatAxisNumber(value) {
  const number = Number(value);
  if (!Number.isFinite(number)) return "?";
  if (Math.abs(number) >= 1000 || (number !== 0 && Math.abs(number) < 0.01)) {
    return number.toExponential(1);
  }
  return String(Math.round(number * 100) / 100);
}''')

# --- 2. the constant -----------------------------------------------------------
sub('''const GENE_SET_MODES = new Set(["dotplot", "combplot"]);''',
    '''const GENE_SET_MODES = new Set(["dotplot", "combplot"]);

// The Coordinates entry that means "not an embedding: use two obs columns".
const OBS_AXES_KEY = "__obs__";''')

# --- 3. hide the axis lists whenever the plot type changes ---------------------
sub('''    if (colorByField) colorByField.classList.toggle("hidden", !wantsUmapOptions);
    if (coordsField) coordsField.classList.toggle("hidden", !wantsUmapOptions);
    if (wantsUmapOptions) refreshUmapOptions(panelKey);''',
    '''    if (colorByField) colorByField.classList.toggle("hidden", !wantsUmapOptions);
    if (coordsField) coordsField.classList.toggle("hidden", !wantsUmapOptions);
    if (wantsUmapOptions) refreshUmapOptions(panelKey);
    syncUmapAxisFields(panelKey);''')

# --- 4. the request -----------------------------------------------------------
sub('''      if (mode === "cluster") {
        const colorBy = getPanelSelectValue(panelKey, "colorby");
        const coordsKey = getPanelSelectValue(panelKey, "coords");
        if (colorBy) params.set("color_by", colorBy);
        if (coordsKey) params.set("coords", coordsKey);
      }''',
    '''      if (mode === "cluster") {
        const colorBy = getPanelSelectValue(panelKey, "colorby");
        const coordsKey = getPanelSelectValue(panelKey, "coords");
        if (colorBy) params.set("color_by", colorBy);
        if (coordsKey === OBS_AXES_KEY) {
          // Two obs columns replace the embedding, so no coords key is sent.
          const xField = getPanelSelectValue(panelKey, "xfield");
          const yField = getPanelSelectValue(panelKey, "yfield");
          if (xField && yField) {
            params.set("x_field", xField);
            params.set("y_field", yField);
          }
        } else if (coordsKey) {
          params.set("coords", coordsKey);
        }
      }''')

# --- 5. a metadata plot needs real axes, and the caption needs the count -------
sub('''  Object.assign(layout, buildSquareUmapAxes(umapData.query, 0.06));
  Plotly.newPlot(panelPlotId(panelKey), traces, layout);
}

function renderPanelExpression(panelKey, expressionData, mode, dotScale) {''',
    '''  if (umapData.axes_source === "obs") {
    // A UMAP hides its tick labels because the numbers mean nothing. A metadata
    // axis is the opposite: the value is the point of the plot, so it is
    // labelled and its ticks are shown.
    layout.xaxis = { title: { text: umapData.x_label, font: { size: 12 } }, showgrid: false, zeroline: false };
    layout.yaxis = { title: { text: umapData.y_label, font: { size: 12 } }, showgrid: false, zeroline: false };
    layout.margin = { t: 20, l: 60, r: 20, b: 55 };
  } else {
    Object.assign(layout, buildSquareUmapAxes(umapData.query, 0.06));
  }
  Plotly.newPlot(panelPlotId(panelKey), traces, layout);
}

function renderPanelExpression(panelKey, expressionData, mode, dotScale) {''')

sub('''  if (umapData && umapData.reference_hidden) {
    parts.push(`Colored by ${umapData.color_label}`);
    parts.push(`Coordinates: ${umapData.coords_label}`);
    parts.push("Reference atlas hidden: it carries neither this annotation nor these coordinates");
  }''',
    '''  if (umapData && umapData.reference_hidden) {
    parts.push(`Colored by ${umapData.color_label}`);
    parts.push(umapData.axes_source === "obs"
      ? `Axes: ${umapData.x_label} against ${umapData.y_label}`
      : `Coordinates: ${umapData.coords_label}`);
    parts.push("Reference atlas hidden: it carries neither this annotation nor these coordinates");
  }
  // Cells the panel could not place. Silence here would read as an absence of
  // cells rather than an absence of a recorded value.
  if (umapData && Number(umapData.n_dropped_no_coordinate) > 0) {
    parts.push(`${umapData.n_dropped_no_coordinate} of ${umapData.n_cells_selected} cells not drawn: no value on both axes`);
  }''')

open(P, "w").write(src)
print("patch14 (app.js axes) applied")
