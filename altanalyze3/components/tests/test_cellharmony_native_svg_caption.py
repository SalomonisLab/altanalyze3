"""The native vector-PDF export keeps its complete legend below the fitted graph."""
from pathlib import Path
import subprocess
import pytest

@pytest.mark.parametrize('landscape', [False, True])
def test_native_svg_pdf_graph_and_bottom_legend_do_not_overlap(landscape):
    source=(Path(__file__).parents[1]/'cellHarmony/webapp/static/app.js').read_text()
    start=source.index('async function saveSvgMarkupAsPdf(')
    end=source.index('\nasync function exportPlotlyElementToPdf(',start)
    function=source[start:end]
    script='''
const assert=require('assert');
const landscape=LANDSCAPE;
const pageWidth=landscape?842:595,pageHeight=landscape?595:842;
let graph,legend;
class Pdf {
 constructor(){this.internal={pageSize:{getWidth:()=>pageWidth,getHeight:()=>pageHeight}};}
 setFontSize(){}
 splitTextToSize(){return Array(10).fill('Complete scientific legend');}
 text(lines,x,y){legend={top:y-9,bottom:y+(lines.length-1)*11};}
 async svg(svg,box){graph=box;}
 save(){}
}
const attrs={width:'1000',height:landscape?'700':'1000'};
const svg={nodeName:'svg',setAttribute:(k,v)=>attrs[k]=v,getAttribute:k=>attrs[k],hasAttribute:k=>k in attrs};
class DOMParser{parseFromString(){return {documentElement:svg};}}
const ensureJsPdfLoaded=async()=>Pdf,ensureSvg2PdfLoaded=async()=>{};
const getSvgIntrinsicSize=()=>({width:Number(attrs.width),height:Number(attrs.height)});
FUNCTION
(async()=>{
 await saveSvgMarkupAsPdf('<svg/>','test.pdf','Full caption');
 // Binary SVG fitting can differ from the exact page margin by ~3e-14 points.
 const geometryRoundoff=1e-12;
 assert(graph.y>=18-geometryRoundoff,JSON.stringify(graph));
 assert(graph.y+graph.height<legend.top,'Graph overlaps or follows its legend');
 assert(legend.bottom<pageHeight-18,'Legend clipped below the page');
 assert(graph.x>=18-geometryRoundoff && graph.x+graph.width<=pageWidth-18+geometryRoundoff,JSON.stringify(graph));
 assert('viewBox' in attrs,'Native scale must retain original SVG coordinates');
})().catch(e=>{console.error(e);process.exit(1);});
'''.replace('LANDSCAPE',str(landscape).lower()).replace('FUNCTION',function)
    result=subprocess.run(['node','-e',script],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
