"""Convierte instancias Minibex en "cintas" (DAG de operaciones, con subexpresiones compartidas y constantes
resueltas) y las cajas de las mediciones, para el evaluador en C++ (aa2.cpp).

Formato (texto, una linea por registro):
  I <instancia> <nvars> <ncajas>
  B <caminata> <prof> lo_1 hi_1 ... lo_n hi_n        (ncajas lineas)
  E <expr> <nnodos>
  <op> <args...>                                       (nnodos lineas; el ultimo nodo es la raiz)
Operaciones: V i | C c | ADD a b | SUB a b | MUL a b | SCL a c | NEG a | REC a | POW a k | EXP a | LOG a |
             SQRT a | SIN a | COS a | ABS a       (a, b: indices de nodos anteriores; MUL a a = cuadrado)
Las cajas son las de medir_cajas.walks con la semilla crc32 (mismas que cajas_m3.csv).
Uso: python3 emitir_cintas.py <salida.txt> <bch>...
"""
import ast, math, sys, zlib
from pathlib import Path
H_ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(H_ / "raak-errores-simbolos"))
import aa_sym as A, medir_cajas as MC  # noqa: E402


class Tape:
    def __init__(self, ev, var_idx):
        self.ev, self.vi, self.nodes, self.memo = ev, var_idx, [], {}

    def add(self, *rec):
        self.nodes.append(rec); return len(self.nodes) - 1

    def emit(self, n):
        key = ast.dump(n)
        if key in self.memo: return self.memo[key]
        ev = self.ev
        if ev.is_const(n): r = self.add("C", repr(ev.const(n)))
        elif isinstance(n, ast.Name): r = self.add("V", self.vi[n.id])
        elif isinstance(n, ast.UnaryOp):
            a = self.emit(n.operand); r = self.add("NEG", a) if isinstance(n.op, ast.USub) else a
        elif isinstance(n, ast.BinOp):
            if isinstance(n.op, (ast.Add, ast.Sub)):
                r = self.add("ADD" if isinstance(n.op, ast.Add) else "SUB", self.emit(n.left), self.emit(n.right))
            elif isinstance(n.op, ast.Mult):
                if ev.is_const(n.left): r = self.add("SCL", self.emit(n.right), repr(ev.const(n.left)))
                elif ev.is_const(n.right): r = self.add("SCL", self.emit(n.left), repr(ev.const(n.right)))
                else: r = self.add("MUL", self.emit(n.left), self.emit(n.right))
            elif isinstance(n.op, ast.Div):
                if ev.is_const(n.right): r = self.add("SCL", self.emit(n.left), repr(1.0 / ev.const(n.right)))
                else: r = self.add("MUL", self.emit(n.left), self.add("REC", self.emit(n.right)))
            elif isinstance(n.op, ast.Pow):
                k = ev.const(n.right)
                if k == int(k) and k > 0:
                    a = self.emit(n.left)
                    r = self.add("MUL", a, a) if k == 2 else (a if k == 1 else self.add("POW", a, int(k)))
                else: r = self.add("EXP", self.add("SCL", self.add("LOG", self.emit(n.left)), repr(k)))
            else: raise ValueError(ast.dump(n))
        elif isinstance(n, ast.Call):
            name = {"exp": "EXP", "log": "LOG", "sqrt": "SQRT", "sin": "SIN", "cos": "COS", "abs": "ABS"}[n.func.id]
            r = self.add(name, self.emit(n.args[0]))
        else: raise ValueError(ast.dump(n))
        self.memo[key] = r
        return r


def main():
    out, files = sys.argv[1], sys.argv[2:]
    with open(out, "w") as fo:
        for f in files:
            box0, consts, trees = A.read_bch(f)
            names = list(box0)
            vi = {v: i for i, v in enumerate(names)}
            boxes = list(MC.walks(box0, seed=zlib.crc32(Path(f).stem.encode())))
            fo.write(f"I {Path(f).stem} {len(names)} {len(boxes)}\n")
            for w, d, b in boxes:
                fo.write(f"B {w} {d} " + " ".join(f"{b[v][0]!r} {b[v][1]!r}" for v in names) + "\n")
            ev = A.Evaluator(box0, consts, "M0")
            for name, tree in trees:
                t = Tape(ev, vi)
                t.emit(tree)
                fo.write(f"E {name} {len(t.nodes)}\n")
                for rec in t.nodes: fo.write(" ".join(str(x) for x in rec) + "\n")


if __name__ == "__main__":
    main()
