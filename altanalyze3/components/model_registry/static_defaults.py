"""Read literal model selectors without importing or executing model code."""
import ast
import posixpath


def api_selection(code, api_path, read_source=None, default_name="DEFAULT_BUNDLE_PATH", depth=0):
    """Interpret only static Path/default constants; never execute historical code."""
    def evaluate(node, constants):
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            return api_path if node.id == "__file__" else constants[node.id]
        if isinstance(node, ast.Dict):
            return {evaluate(k, constants): evaluate(v, constants) for k, v in zip(node.keys, node.values)}
        if isinstance(node, ast.Subscript):
            return evaluate(node.value, constants)[evaluate(node.slice, constants)]
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
            return posixpath.join(evaluate(node.left, constants), evaluate(node.right, constants))
        if isinstance(node, ast.Attribute) and node.attr == "parent":
            return posixpath.dirname(evaluate(node.value, constants))
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name) and node.func.id == "Path":
                return evaluate(node.args[0], constants)
            if isinstance(node.func, ast.Attribute) and node.func.attr == "resolve":
                return posixpath.normpath(evaluate(node.func.value, constants))
            if isinstance(node.func, ast.Attribute) and node.func.attr == "with_name":
                return posixpath.join(posixpath.dirname(evaluate(node.func.value, constants)), evaluate(node.args[0], constants))
        raise ValueError("Not a static default expression")
    constants = {}
    for node in ast.parse(code).body:
        if isinstance(node, ast.ImportFrom) and node.level == 1 and node.module and read_source and depth < 8:
            imported_path = posixpath.join(posixpath.dirname(api_path), node.module.replace(".", "/") + ".py")
            imported_code = read_source(imported_path)
            if imported_code:
                for alias in node.names:
                    value = api_selection(imported_code, imported_path, read_source, alias.name, depth + 1)
                    if value is not None:
                        constants[alias.asname or alias.name] = value
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                if isinstance(target, ast.Name):
                    try:
                        constants[target.id] = evaluate(node.value, constants)
                    except (ValueError, KeyError, TypeError, AttributeError):
                        pass
    return constants.get(default_name)

