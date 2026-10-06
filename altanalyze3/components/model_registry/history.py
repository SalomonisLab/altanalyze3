"""Export committed scALABLE default selections and artifact Git identities.

Commit times describe repository states, not deployment dates or training times.
Only Git objects are read; historical models are never loaded or inferred from names.
"""
import json
import posixpath
import subprocess
from pathlib import Path, PurePosixPath
from .registry import PACKAGE_ROOT
from .static_defaults import api_selection


REPO_ROOT = subprocess.check_output(["git", "-C", str(PACKAGE_ROOT), "rev-parse", "--show-toplevel"], text=True).strip()


def git(*args):
    return subprocess.check_output(["git", "-C", REPO_ROOT, *args], text=True).strip()




def main():
    prefix = PACKAGE_ROOT.relative_to(REPO_ROOT).as_posix() + "/"
    config_path = prefix + "components/cellHarmony/flask/reference_config.json"
    tracked_paths = [config_path] + [prefix + "components/" + name for name in
        ("rna2lipid", "rna2adt", "rna2grn", "rna2metabolite", "fastComm")]
    commits = git("log", "--reverse", "--format=%H %cI", "--", *tracked_paths).splitlines()
    history = []
    previous = None
    for line in commits:
        commit, timestamp = line.split(" ", 1)
        objects = {parts[1]: parts[0].split()[2] for parts in
                   (line.split("\t", 1) for line in git("ls-tree", "-r", commit).splitlines())}
        config = json.loads(git("show", f"{commit}:{config_path}")) if config_path in objects else {}
        selections = []
        for species in config.get("species", []):
            for ref in species.get("references", []):
                for modality, cfg in ref.get("impute_config", {}).items():
                    if "bundle_path" not in cfg:
                        continue
                    import posixpath
                    path = posixpath.normpath(str(PurePosixPath(config_path).parent / cfg["bundle_path"]))
                    selections.append({"reference": ref["id"], "modality": modality,
                        "config": cfg, "artifact_git_blob": objects.get(path), "artifact_path": path})
        for component in ("rna2lipid", "rna2adt", "rna2grn", "rna2metabolite", "rna2lipid/aml"):
            api_path = prefix + "components/" + component + "/api.py"
            if api_path not in objects:
                continue
            path = api_selection(git("show", f"{commit}:{api_path}"), api_path,
                                 lambda path: git("show", f"{commit}:{path}") if path in objects else None)
            if path:
                selections.append({"application": "altanalyze3", "context": "api-default",
                    "component": component, "artifact_path": path,
                    "artifact_git_blob": objects.get(path), "api_git_blob": objects[api_path]})
        if selections != previous:
            if history:
                history[-1]["repository_until"] = timestamp
            history.append({"commit": commit, "repository_from": timestamp,
                            "repository_until": None, "defaults": selections})
            previous = selections
    destination = PACKAGE_ROOT / "model_repository/default_history.json"
    destination.write_text(json.dumps({"schema_version": "1.0",
        "scope": "committed AltAnalyze3 API and scALABLE reference selections; not deployment history",
        "history": history}, indent=2, sort_keys=True) + "\n")
    print(f"Exported {len(history)} committed default intervals")


if __name__ == "__main__":
    main()
