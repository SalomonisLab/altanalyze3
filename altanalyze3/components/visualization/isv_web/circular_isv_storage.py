"""Disk-backed caches for the unchanged circRNA ISV display selection.

Collections preserve source insertion order. They are implementation caches;
source evidence and selected tasks remain authoritative and unchanged.
"""
import json, sqlite3
from itertools import groupby

def packed(value):
    return json.dumps(value,ensure_ascii=False,separators=(',',':'))

def unpacked(value):
    value=json.loads(value)
    return tuple(value) if isinstance(value,list) else value

class Store:
    def __init__(self,path):
        self.db=sqlite3.connect(path)
        self.db.execute('PRAGMA journal_mode=OFF')
        self.db.execute('PRAGMA synchronous=OFF')
        self.db.execute('PRAGMA temp_store=FILE')
        self.db.execute('PRAGMA cache_size=-32768')
    def commit(self):self.db.commit()
    def close(self):self.db.commit();self.db.close()

class DiskMap:
    def __init__(self,store,name):
        assert name.isidentifier();self.db=store.db;self.name=name
        self.db.execute(f'CREATE TABLE {name}(key TEXT PRIMARY KEY,value TEXT)')
    def __setitem__(self,key,value):
        self.db.execute(f'INSERT OR REPLACE INTO {self.name} VALUES (?,?)',(packed(key),packed(value)))
    def get(self,key,default=None):
        r=self.db.execute(f'SELECT value FROM {self.name} WHERE key=?',(packed(key),)).fetchone()
        return json.loads(r[0]) if r else default
    def __getitem__(self,key):
        r=self.db.execute(f'SELECT value FROM {self.name} WHERE key=?',(packed(key),)).fetchone()
        if r is None:raise KeyError(key)
        return json.loads(r[0])

class DiskSet:
    def __init__(self,store,name):
        assert name.isidentifier();self.db=store.db;self.name=name
        self.db.execute(f'CREATE TABLE {name}(key TEXT PRIMARY KEY)')
    def add(self,key):self.db.execute(f'INSERT OR IGNORE INTO {self.name} VALUES (?)',(packed(key),))
    def __contains__(self,key):return self.db.execute(f'SELECT 1 FROM {self.name} WHERE key=?',(packed(key),)).fetchone() is not None

class DiskSetMap:
    def __init__(self,store,name):
        assert name.isidentifier();self.db=store.db;self.name=name
        self.db.execute(f'CREATE TABLE {name}(key TEXT,value TEXT,PRIMARY KEY(key,value))')
    def add(self,key,value):self.db.execute(f'INSERT OR IGNORE INTO {self.name} VALUES (?,?)',(packed(key),packed(value)))
    def get(self,key,default=None):
        found={json.loads(r[0]) for r in self.db.execute(f'SELECT value FROM {self.name} WHERE key=?',(packed(key),))}
        return found if found else default
    def __getitem__(self,key):
        r=self.get(key)
        if r is None:raise KeyError(key)
        return r
    def items(self):
        for key,group in groupby(self.db.execute(f'SELECT key,value FROM {self.name} ORDER BY key'),lambda r:r[0]):
            yield unpacked(key),{json.loads(r[1]) for r in group}

class DiskPairs:
    def __init__(self,store):
        self.db=store.db
        self.db.execute('CREATE TABLE pair_groups(ordinal INTEGER PRIMARY KEY,event TEXT,library TEXT,query TEXT,UNIQUE(event,library,query))')
        self.db.execute('CREATE TABLE pairs(ordinal INTEGER PRIMARY KEY,group_id INTEGER,payload TEXT)')
        self.db.execute('CREATE INDEX pair_group ON pairs(group_id,ordinal)')
        self.db.execute('CREATE INDEX pair_event ON pair_groups(event,ordinal)')
    def add(self,key,value):
        self.db.execute('INSERT OR IGNORE INTO pair_groups(event,library,query) VALUES (?,?,?)',key)
        gid=self.db.execute('SELECT ordinal FROM pair_groups WHERE event=? AND library=? AND query=?',key).fetchone()[0]
        self.db.execute('INSERT INTO pairs(group_id,payload) VALUES (?,?)',(gid,packed(value)))
    def items(self,event=None):
        sql='SELECT ordinal,event,library,query FROM pair_groups'
        args=()
        if event is not None:sql+=' WHERE event=?';args=(event,)
        for gid,ev,lib,q in self.db.execute(sql+' ORDER BY ordinal',args):
            yield (ev,lib,q),[json.loads(r[0]) for r in self.db.execute('SELECT payload FROM pairs WHERE group_id=? ORDER BY ordinal',(gid,))]

def write_json_array(path,values):
    with open(path,'w') as h:
        h.write('[');sep=''
        for value in values:h.write(sep);json.dump(value,h);sep=','
        h.write(']')
