P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/templates/index.html"
src = open(P).read()

for panel in ("viz1", "viz2"):
    old = f'''                      <label class="field hidden inline-field compact-inline-field" id="{panel}-coords-field">
                        <span>Coordinates</span>
                        <select id="{panel}-coords" class="compact-select"></select>
                      </label>'''
    new = old + f'''
                      <!-- Axes. Choosing "obs columns" in Coordinates puts any two
                           numeric cell annotations on the axes, the way ShinyCell
                           plots one metadata field against another. Counts live in
                           the matrix, not in obs, so nothing here is expression. -->
                      <label class="field hidden inline-field compact-inline-field" id="{panel}-xfield-field">
                        <span>X axis</span>
                        <select id="{panel}-xfield" class="compact-select"></select>
                      </label>
                      <label class="field hidden inline-field compact-inline-field" id="{panel}-yfield-field">
                        <span>Y axis</span>
                        <select id="{panel}-yfield" class="compact-select"></select>
                      </label>'''
    assert src.count(old) == 1, f"{panel}: {src.count(old)}"
    src = src.replace(old, new)

open(P, "w").write(src)
print("patch13 (index.html axes) applied")
