# Python 2.7 / 3 compatible adapter for the unmodified, externally installed AltAnalyze2.
import os, sys, runpy, json, collections
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.fonttype']='none'
import matplotlib.pyplot as pylab

def main(config_path):
    cfg=json.load(open(config_path));base=cfg['base'];expected=collections.Counter()
    for line in open(base+'/coordinates.tsv'):
        t=line.rstrip().split('\t');expected[(float(t[1]),float(t[2]))]+=1
    original_save=pylab.savefig
    def checked_save(path,*args,**kwargs):
        if str(path).endswith('.pdf'):
            seen=collections.Counter();groups={}
            for ax in pylab.gcf().axes:
                for coll in ax.collections:
                    if hasattr(coll,'get_offsets'):
                        xy=coll.get_offsets()
                        if len(xy):
                            groups[coll.get_label()]=len(xy)
                            seen.update((float(v[0]),float(v[1])) for v in xy)
            if seen!=expected:raise ValueError('Native plotted coordinates differ from the complete saved UMAP')
            if groups!=cfg['expected_groups']:raise ValueError('Native group counts differ')
            json.dump(dict(native_point_count=sum(seen.values()),exact_coordinate_multiset_equal=True,native_group_sizes=groups),open(base+'/native_plot_validation.json','w'),indent=2)
            original_save(str(path)[:-4]+'.svg')
        return original_save(path,*args,**kwargs)
    pylab.savefig=checked_save
    app=cfg['app']
    if sys.version_info[0] < 3 and isinstance(app,unicode):app=app.encode('utf-8')
    sys.path.insert(0,os.path.dirname(app))
    sys.argv=[app,'--image','UMAP','--plotType','2D','--display','False','--species','Hs','--platform','RNASeq','--zscore','yes','--input',base+'/exp.KINNEX-clones.txt','--coordinateFile',base+'/coordinates.tsv','--reimportModelScores','True']
    if sys.version_info[0] < 3:sys.argv=[v.encode('utf-8') if isinstance(v,unicode) else v for v in sys.argv]
    try:runpy.run_path(app,run_name='__main__')
    finally:pylab.savefig=original_save
    if not os.path.isfile(base+'/native_plot_validation.json'):raise ValueError('Native plot not verified')
if __name__=='__main__':main(sys.argv[1])
