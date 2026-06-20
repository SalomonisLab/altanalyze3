"""Parse a WikiPathways GPML release (zip) into pathway gene sets (Ensembl ENSG) and curated
interaction topology, for pathway-constrained network contrasts."""

import re
import zipfile
import logging
import xml.etree.ElementTree as ET

from .. import config

logger = logging.getLogger(__name__)


def _localname(tag):
    return tag.rsplit("}", 1)[-1]


def load_pathways(zip_path=None):
    """Return {wpid: {'name':..., 'genes': set(ENSG), 'interactions': [(srcENSG, dstENSG), ...]}}.
    DataNodes are joined by their Ensembl Xref (genes keyed by ENSG, matching our ENSG column);
    complexes are expanded to member genes via GraphId; metabolite/other nodes are skipped."""
    zip_path = zip_path or config.WIKIPATHWAYS_ZIP
    pathways = {}
    with zipfile.ZipFile(zip_path) as z:
        for member in z.namelist():
            if not member.endswith(".gpml"):
                continue
            try:
                root = ET.fromstring(z.read(member))
            except ET.ParseError:
                continue
            name = root.attrib.get("Name", "")
            m = re.search(r"(WP\d+)", member)
            wpid = m.group(1) if m else member
            graphid_to_ensg = {}
            genes = set()
            for node in root.iter():
                if _localname(node.tag) != "DataNode":
                    continue
                gid = node.attrib.get("GraphId", "")
                ensg = None
                for child in node:
                    if _localname(child.tag) == "Xref" and child.attrib.get("Database") == "Ensembl":
                        xid = child.attrib.get("ID", "")
                        if xid.startswith("ENSG"):
                            ensg = xid.split(".")[0]
                if ensg:
                    genes.add(ensg)
                    if gid:
                        graphid_to_ensg[gid] = ensg
            interactions = []
            for inter in root.iter():
                if _localname(inter.tag) != "Interaction":
                    continue
                pts = [p for p in inter.iter() if _localname(p.tag) == "Point"]
                refs = [p.attrib.get("GraphRef", "") for p in pts]
                if len(refs) >= 2 and refs[0] in graphid_to_ensg and refs[-1] in graphid_to_ensg:
                    interactions.append((graphid_to_ensg[refs[0]], graphid_to_ensg[refs[-1]]))
            if genes:
                pathways[wpid] = {"name": name, "genes": genes, "interactions": interactions}
    logger.info("[wikipathways] parsed %d pathways (with >=1 Ensembl gene)", len(pathways))
    return pathways
