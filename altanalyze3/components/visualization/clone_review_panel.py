"""Assemble saved SVG artwork on a half-height US Letter artboard with >=5 pt editable text.

Only figure layout changes. Plot paths/images remain sourced from the supplied
standalone figures; labels are transformed independently to avoid font shrinkage.
"""
from pathlib import Path
import copy,json,re,argparse,xml.etree.ElementTree as ET
from playwright.sync_api import sync_playwright
NS='http://www.w3.org/2000/svg';ET.register_namespace('',NS)
XL='http://www.w3.org/1999/xlink';ET.register_namespace('xlink',XL)

def source_panel(file,groups,prefix):
 r=ET.parse(file).getroot();fig=r.find('.//*[@id="figure_1"]')
 for e in list(fig):
  if e.attrib.get('id') not in groups:fig.remove(e)
 if prefix=='c':
  for g in r.findall('.//{'+NS+'}g'):
   if len(g.findall('{'+NS+'}path'))==256:g.set('shape-rendering','crispEdges')
 ids={e.attrib['id']:prefix+'_'+e.attrib['id'] for e in r.iter() if 'id' in e.attrib}
 for e in r.iter():
  for k,v in list(e.attrib.items()):
   if k=='id':e.set(k,ids[v])
   else:
    v=re.sub(r'url\(#([^)]*)\)',lambda m:'url(#'+ids.get(m[1],m[1])+')',v)
    if v.startswith('#') and v[1:] in ids:v='#'+ids[v[1:]]
    e.set(k,v)
 return ET.tostring(r,encoding='unicode')

JS=r'''({prefix,w,h,uniform})=>{
 const root=document.querySelector('svg'),group=root.querySelector('[id="'+prefix+'_figure_1"]');
 const box=group.getBBox();let sx=w/box.width,sy=h/box.height;
 if(uniform)sx=sy=Math.min(sx,sy);
 const original=group.cloneNode(true);
 const arrowRecords=prefix==='b'?[...group.querySelectorAll(':scope > g > g[id*="patch_"]')].filter(g=>g.querySelectorAll(':scope > path').length===2).map(g=>{const v=g.querySelector('path').getAttribute('d').match(/[-+]?[0-9]*\.?[0-9]+(?:e[-+]?[0-9]+)?/gi).map(Number);return {id:g.id,x1:v[0],y1:v[1],x2:v[v.length-2],y2:v[v.length-1]};}):[];
 const circleBoxes=prefix==='b'?[...group.querySelectorAll('[id*="AnnotationBbox"]')].map(g=>{const b=g.getBBox();return {id:g.id,cx:b.x+b.width/2,cy:b.y+b.height/2};}):[];
 let records=[...group.querySelectorAll('text')].map(t=>{
  const m=root.getCTM().inverse().multiply(t.getCTM());
  const pt=new DOMPoint(parseFloat(t.getAttribute('x')||0),parseFloat(t.getAttribute('y')||0)).matrixTransform(m);
  const st=getComputedStyle(t);const gb=t.parentNode.getBBox();const gc=new DOMPoint(gb.x+gb.width/2,gb.y).matrixTransform(root.getCTM().inverse().multiply(t.parentNode.getCTM()));
  return {gid:t.parentNode.id,gx:gc.x,text:t.textContent,x:pt.x,y:pt.y,size:parseFloat(st.fontSize),style:t.getAttribute('style')||'',anchor:st.textAnchor,
   angle:Math.atan2(m.b,m.a)*180/Math.PI};
 });
 if(prefix==='b'){
  const ids=[...new Set(records.map(r=>r.gid))];
  for(const id of ids){
   const rows=records.filter(r=>r.gid===id);
   if(rows.length && rows[0].text.startsWith('+') && rows.map(r=>r.text).join('').includes('c.')){
    const first=rows[0],joined=rows.map(r=>r.text).join(''),split=joined.indexOf('c.'),tokens=[joined.slice(0,split).trim(),...joined.slice(split).replace(/\s+/g,'').match(/.{1,12}/g)];
    records=records.filter(r=>r.gid!==id);
    tokens.forEach((v,i)=>records.push({...first,text:v,x:first.gx,y:first.y+i*first.size*1.1,anchor:'middle'}));
   }
  }
 }
 if(prefix==='c'||prefix==='d'){
  records=records.filter(r=>!(prefix==='c' && (/^WT [0-9]/.test(r.text)||r.text==='(schematic)')));
  for(const r of records)r.text=({'Prediagnosis':'Initiation','Diagnosis / pre-HMA':'Diagnosis','sAML relapse':'Relapse','sAML CR':'sAML-CR'})[r.text]||r.text;
 }
 let final,b,fonts=[],omitted=[];
 for(let iteration=0;iteration<8;iteration++){
  final=document.createElementNS(root.namespaceURI,'g');final.id='layout';
  const art=original.cloneNode(true);art.querySelectorAll('text').forEach(t=>t.remove());
  if(prefix==='b')for(const q of circleBoxes){const g=art.querySelector('[id="'+q.id+'"]');const z=Math.min(sx,sy);g.setAttribute('transform',`translate(${q.cx} ${q.cy}) scale(${z/sx} ${z/sy}) translate(${-q.cx} ${-q.cy})`);}
  if(prefix==='b'){
   for(const q of arrowRecords)art.querySelector('[id="'+q.id+'"]').remove();
   for(const path of art.querySelectorAll('path'))if(/fill:\s*#ffffff(?:;|$)/.test(path.getAttribute('style')||''))path.remove();
  }
  if(prefix==='d')for(const path of art.querySelectorAll('path'))if(/stroke:\s*#000000/.test(path.getAttribute('style')||'')){path.style.strokeWidth='0.2';path.setAttribute('vector-effect','non-scaling-stroke');}
  art.setAttribute('transform',`translate(${-box.x*sx} ${-box.y*sy}) scale(${sx} ${sy})`);final.append(art);fonts=[];
  if(prefix==='b')for(const q of arrowRecords){
   const x1=(q.x1-box.x)*sx,y1=(q.y1-box.y)*sy,x2=(q.x2-box.x)*sx,y2=(q.y2-box.y)*sy;
   const theta=Math.atan2(y2-y1,x2-x1),len=2.5,spread=Math.PI/7;
   const path=document.createElementNS(root.namespaceURI,'path');
   path.setAttribute('d',`M ${x1} ${y1} L ${x2} ${y2} M ${x2-len*Math.cos(theta-spread)} ${y2-len*Math.sin(theta-spread)} L ${x2} ${y2} L ${x2-len*Math.cos(theta+spread)} ${y2-len*Math.sin(theta+spread)}`);
   path.setAttribute('style','fill:none;stroke:#444444;stroke-width:0.7;stroke-linecap:round;stroke-linejoin:round');path.setAttribute('data-role','tree-arrow');final.append(path);
  }
  const labelY=new Map();
  if(prefix==='b')for(const id of [...new Set(records.map(r=>r.gid))]){
   const rows=records.filter(r=>r.gid===id).sort((a,b)=>a.y-b.y);
   rows.forEach((r,i)=>labelY.set(r,(rows[0].y-box.y)*sy+i*Math.max(5.6,rows.length>1?(rows[1].y-rows[0].y)*sy:0)));
  }
  if(prefix==='d'){
   const rows=records.filter(r=>r.anchor==='end' && !/^n\s*=|^[0-9]+$/.test(r.text)).sort((a,b)=>a.y-b.y);
   let next=h-5;
   for(const r of [...rows].reverse()){const y=Math.min((r.y-box.y)*sy,next);labelY.set(r,y);next=y-5.7;}
  }
  for(const r of records){
   const t=document.createElementNS(root.namespaceURI,'text');let x=(r.x-box.x)*sx,y=labelY.has(r)?labelY.get(r):(r.y-box.y)*sy;if(prefix==='c' && Math.abs(r.angle)>80)x-=7;
   const stage=prefix==='d'&&/^(Pre-HMA|On-HMA|Post-HMA|Diagnosis|sAML|sAML-CR|Relapse)$/.test(r.text);const size=stage?5.05:Math.max(5.05,Math.min(prefix==='b'?5.5:Infinity,r.size*Math.min(sx,sy)));fonts.push(size);
   t.textContent=r.text;t.setAttribute('style',r.style);t.style.fontSize=size+'px';t.style.fontFamily='Arial';if(stage)t.style.fontWeight='normal';
   t.setAttribute('x',x);t.setAttribute('y',y);t.style.textAnchor=r.anchor;
   if(Math.abs(r.angle)>0.01)t.setAttribute('transform',`rotate(${r.angle} ${x} ${y})`);
   final.append(t);
  }
  const old=root.querySelector('#layout');if(old)old.remove();root.append(final);group.setAttribute('display','none');
  const kept=[];omitted=[];
  for(const t of final.querySelectorAll(':scope > text')){
   if(!/^(SC[0-9]+|WT)$/.test(t.textContent))continue;
   const q=t.getBBox();
   if(kept.some(b=>q.x<b.x+b.width+0.7&&q.x+q.width+0.7>b.x&&q.y<b.y+b.height&&q.y+q.height>b.y)){
    omitted.push(t.textContent);t.remove();
   }else kept.push(q);
  }
  b=final.getBBox();if(b.width<=w+.01&&b.height<=h+.01)break;
  let fx=Math.min(1,w/b.width*.995),fy=Math.min(1,h/b.height*.995);
  if(uniform)fx=fy=Math.min(fx,fy);sx*=fx;sy*=fy;
 }
 final.setAttribute('transform',`translate(${(w-b.width)/2-b.x} ${(h-b.height)/2-b.y})`);
 return {xml:final.outerHTML,defs:[...root.querySelectorAll(':scope > defs')].map(d=>d.outerHTML).join(''),width:b.width,height:b.height,fonts,sx,sy,omitted};
}'''

def compose(figures,output,patient=None,chrome=None):
 if chrome is None:
  import os,shutil
  chrome=os.environ.get('ALTANALYZE_CHROME') or shutil.which('chromium') or shutil.which('google-chrome')
  if not chrome and Path('/Applications/Google Chrome.app/Contents/MacOS/Google Chrome').exists():chrome='/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
 F=Path(figures);O=Path(output);O.mkdir(parents=True,exist_ok=True);patient=patient or F.parent.name
 include_unassigned=json.loads((F/'display_reporting.json').read_text()).get('include_unassigned',False) if (F/'display_reporting.json').exists() else False
 heat=F/'heatmap'/f'{patient}_compact_source_heatmap_thicker.svg'
 if not heat.exists():heat=F/'compact_source_heatmap'/f'{patient}_compact_source_heatmap_thicker.svg'
 alluvial=F/'alluvial/broad_state_clonal_alluvial.svg'
 if not alluvial.exists():alluvial=next((F/'state_clone_alluvial').glob('*alluvial.svg'))
 if chrome is None and Path('/Applications/Google Chrome.app/Contents/MacOS/Google Chrome').exists():chrome='/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
 t=ET.parse(heat);axs=[g.attrib['id'] for g in t.findall('.//{'+NS+'}g') if re.fullmatch('axes_[0-9]+',g.attrib.get('id',''))]
 heatgroups=axs[:-1] # last axes is the duplicated full-genotype key, supplied by tree/companion tables
 specs=[('a',heat,heatgroups,(18,28,350,140),False),
        ('b',F/'clonal_tree.svg',['axes_1'],(18,198,350,108),False),
        ('b_key',F/'clonal_tree.svg',['axes_2'],(18,311,350,37),True),
        ('c',F/'fishplot.svg',['axes_1'],(386,28,208,140),False),
        ('d',alluvial,['axes_1'],(386,198,208,150),False)]
 root=ET.Element('{'+NS+'}svg',{'width':'612pt','height':'396pt','viewBox':'0 0 612 396'})
 ET.SubElement(root,'{'+NS+'}rect',{'width':'612','height':'396','fill':'white'})
 audit=[]
 with sync_playwright() as pw:
  browser=pw.chromium.launch(headless=True,**({'executable_path':chrome} if chrome else {}));page=browser.new_page(viewport={'width':816,'height':528},device_scale_factor=2)
  for letter,path,groups,(x,y,w,h),uniform in specs:
   page.set_content(source_panel(path,groups,letter));page.evaluate('document.fonts.ready')
   result=page.evaluate(JS,dict(prefix=letter,w=w,h=h,uniform=uniform))
   if result['width']>w+.1 or result['height']>h+.1:raise ValueError('Panel does not fit: '+letter)
   wrapper=ET.SubElement(root,'{'+NS+'}g',{'id':'panel_'+letter,'transform':f'translate({x} {y})'})
   if result['defs']:
    for defs in ET.fromstring('<svg xmlns="'+NS+'">'+result['defs']+'</svg>'):root.append(defs)
   wrapper.append(ET.fromstring('<svg xmlns="'+NS+'" xmlns:xlink="'+XL+'">'+result['xml']+'</svg>')[0])
   audit.append({'panel':letter,'source':str(path),'placement_pt':[x,y,w,h],'min_font_pt':min(result['fonts']),'text_count':len(result['fonts']),'scale_x':result['sx'],'scale_y':result['sy'],'overlapping_clone_labels_omitted':result['omitted']})
  def text(x,y,txt,size=6,**kw):
   el=ET.SubElement(root,'{'+NS+'}text',{'x':str(x),'y':str(y),'font-family':'Arial','font-size':str(size),**kw});el.text=txt
  text(306,13,'Patient: '+patient,9,**{'text-anchor':'middle'})
  for letter,x,y in [('a',18,24),('b',18,194),('c',386,24),('d',386,194)]:text(x,y,letter,8,**{'font-weight':'bold'})
  # Shared evidence key below panel a: two rows, at publication-size text.
  items=[('Real','#303030'),('PV','#d1352b'),('MPV','#2f6fb7'),('LOH','#2e9a62'),('CNV','#e08a1e'),('DEG','#7b4ea3'),('DEI','#1f9aa6'),('Accepted','#a6a6a6')]
  for i,(lab,col) in enumerate(items):
   x=30+i*41;ET.SubElement(root,'{'+NS+'}rect',{'x':str(x),'y':'173','width':'5','height':'5','fill':col});text(x+7,178,lab,5.5)
  for i,(lab,col) in enumerate([('No call','#ffffff'),('WT/reference','#acd9a5'),('MUT','#f06b48')]):
   x=99+i*72;ET.SubElement(root,'{'+NS+'}rect',{'x':str(x),'y':'183','width':'5','height':'5','fill':col,'stroke':'#aaa','stroke-width':'.4'});text(x+7,188,lab,5.5)
  import pandas as pd,textwrap
  key=pd.read_csv(F/'clone_key.tsv',sep='\t');key=key if include_unassigned else key[key.clone.ne('U')]
  for i,row in enumerate(key.itertuples()):
   x=18+i*52;ET.SubElement(root,'{'+NS+'}rect',{'x':str(x),'y':'353','width':'5','height':'5','fill':row.color});text(x+7,358,row.clone,6)
  cap=('a, Direct and imputed mutation evidence by cell and specimen; white, no call; green, reference; orange, mutant. '
       'b, Saved fitted clonal tree; dashed nodes are latent model states, and arrows do not establish observed mutation-gain order. '
       'c, Minimum supported lineage frequencies across clinical specimens; prediagnosis is schematic. Matched grafts are shown in a. '
       'd, Cell-state/lineage frequencies; left colors identify cell states and right colors identify clones. Connectors interpolate proportions, not tracked cells. '
       'U cells are excluded from displayed panels and denominators; WT (no accepted mutant call) remains. Source assignments and fitted trees are unchanged. '
       'One patient per page; cells are not biological replicates. No hypothesis tests or confidence intervals. Sparse RNA evidence and imputation limit genotype certainty.')
  if include_unassigned:cap=cap.replace('U cells are excluded from displayed panels and denominators; WT (no accepted mutant call) remains.', 'All cells are included; U denotes unresolved lineage assignment and WT denotes no accepted mutant call.')
  cap += ' Small clone labels are omitted where they overlap; all bands remain.'
  for i,line in enumerate(textwrap.wrap(cap,215)):text(18,368+i*6,line,5.05)
  file=O/'review_panel.svg';ET.ElementTree(root).write(file,encoding='utf-8',xml_declaration=True)
  page.set_content('<html><head><style>@page{size:612pt 396pt;margin:0}html,body{margin:0;padding:0;width:612pt;height:396pt}svg{display:block}</style></head><body>'+file.read_text()+'</body></html>');page.evaluate('document.fonts.ready')
  page.screenshot(path=str(O/'review_panel.png'),full_page=False)
  page.pdf(path=str(O/'review_panel.pdf'),width='8.5in',height='5.5in',margin={k:'0' for k in ['top','bottom','left','right']},print_background=True)
  # CTM scales are converted from browser CSS pixels to physical points.
  stats=page.evaluate('''()=>[...document.querySelectorAll('svg text')].map(t=>{const m=t.getScreenCTM(),b=t.getBoundingClientRect();return {text:t.textContent,size_pt:parseFloat(getComputedStyle(t).fontSize)*Math.hypot(m.a,m.b)*.75,x:b.x*.75,y:b.y*.75,w:b.width*.75,h:b.height*.75}})''')
  if min(r['size_pt'] for r in stats)<4.999:raise ValueError('Text below 5 pt')
  if any(r['x']<-.1 or r['y']<-.1 or r['x']+r['w']>612.1 or r['y']+r['h']>396.1 for r in stats):raise ValueError('Text outside page')
  browser.close()
 report={'page_pt':[612,396],'minimum_text_pt':min(r['size_pt'] for r in stats),'panels':audit,'all_text_inside_page':True,'text':stats}
 (O/'review_panel_layout_validation.json').write_text(json.dumps(report,indent=2))
 return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--figures',required=True);p.add_argument('--output',required=True);p.add_argument('--patient');p.add_argument('--chrome');a=p.parse_args();compose(a.figures,a.output,a.patient,a.chrome)
