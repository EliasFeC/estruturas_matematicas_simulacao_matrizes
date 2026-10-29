import math
import random
import time

TOL = 1e-9  # tolerância relativa para comparar números de ponto flutuante


# =====================================================================
# PARTE 0 - Ordens 1 e 2 (fórmulas diretas)
# =====================================================================
def determinante_12(A):
    """Ordem 1: det = a.  Ordem 2: det = a*d - b*c."""
    if len(A) == 1:
        return A[0][0]
    a, b = A[0]
    c, d = A[1]
    return a * d - b * c


def inversa_12(A):
    """Ordem 1: [[1/a]].  Ordem 2: (1/det) * [[d, -b], [-c, a]]. None se det=0."""
    det = determinante_12(A)
    if det == 0:
        return None
    if len(A) == 1:
        return [[1 / det]]
    a, b = A[0]
    c, d = A[1]
    return [[d / det, -b / det], [-c / det, a / det]]

def sarrus_3x3(A, mostrar_passos=True):
    """
    Repete as duas primeiras colunas à direita da matriz:
        a b c | a b
        d e f | d e
        g h i | g h
    Diagonais principais (+): aei, bfg, cdh
    Diagonais secundárias (-): ceg, afh, bdi
    det = (aei + bfg + cdh) - (ceg + afh + bdi)
    """
    a, b, c = A[0]
    d, e, f = A[1]
    g, h, i = A[2]

    pos = [a * e * i, b * f * g, c * d * h]
    neg = [c * e * g, a * f * h, b * d * i]
    det = sum(pos) - sum(neg)

    if mostrar_passos:
        print("\nPasso 1 - repetir as duas primeiras colunas à direita:")
        for lin in A:
            print("   | " + "  ".join(fmt(x).rjust(5) for x in lin)
                  + " |  " + "  ".join(fmt(x).rjust(5) for x in lin[:2]))
        print("\nPasso 2 - diagonais principais:")
        print(f"   {fmt(a)}*{fmt(e)}*{fmt(i)} = {fmt(pos[0])}")
        print(f"   {fmt(b)}*{fmt(f)}*{fmt(g)} = {fmt(pos[1])}")
        print(f"   {fmt(c)}*{fmt(d)}*{fmt(h)} = {fmt(pos[2])}")
        print(f"   soma = {fmt(sum(pos))}")
        print("\nPasso 3 - diagonais secundárias:")
        print(f"   {fmt(c)}*{fmt(e)}*{fmt(g)} = {fmt(neg[0])}")
        print(f"   {fmt(a)}*{fmt(f)}*{fmt(h)} = {fmt(neg[1])}")
        print(f"   {fmt(b)}*{fmt(d)}*{fmt(i)} = {fmt(neg[2])}")
        print(f"   soma = {fmt(sum(neg))}")
        print(f"\nPasso 4 - det(A) = {fmt(sum(pos))} - {fmt(sum(neg))} = {fmt(det)}")
    return det


# =====================================================================
# PARTE 1 - Ordem n qualquer (eliminação de Gauss com pivotamento parcial)
# =====================================================================
def _tolerancia(A):
    """Tolerância de singularidade proporcional à magnitude da matriz."""
    maior = max((abs(x) for linha in A for x in linha), default=1.0)
    return 1e-12 * max(1.0, maior)


def determinante_n(A):
    """
    det por eliminação de Gauss, O(n^3).
    Pivotamento parcial (maior valor da coluna) -> estabilidade numérica.
    Cada troca de linhas inverte o sinal; det = produto da diagonal.
    """
    n = len(A)
    M = [linha[:] for linha in A]
    tol = _tolerancia(A)
    det = 1.0
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(M[i][k]))
        if abs(M[p][k]) < tol:
            return 0.0
        if p != k:
            M[k], M[p] = M[p], M[k]
            det = -det
        det *= M[k][k]
        Mk = M[k]
        for i in range(k + 1, n):
            Mi = M[i]
            f = Mi[k] / Mk[k]
            if f != 0:
                for j in range(k + 1, n):
                    Mi[j] -= f * Mk[j]
    return det


def inversa_n(A):
    """Inversa por Gauss-Jordan sobre [A | I], com pivotamento parcial. None se singular."""
    n = len(A)
    tol = _tolerancia(A)
    M = [list(map(float, A[i])) + [1.0 if i == j else 0.0 for j in range(n)]
         for i in range(n)]
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(M[i][k]))
        if abs(M[p][k]) < tol:
            return None
        M[k], M[p] = M[p], M[k]
        piv = M[k][k]
        M[k] = [x / piv for x in M[k]]
        for i in range(n):
            if i != k:
                f = M[i][k]
                if f != 0:
                    M[i] = [x - f * y for x, y in zip(M[i], M[k])]
    return [linha[n:] for linha in M]


def determinante_laplace(A):
    """Expansão de Laplace recursiva, O(n!) - só para mostrar a limitação."""
    n = len(A)
    if n == 1:
        return A[0][0]
    if n == 2:
        return A[0][0] * A[1][1] - A[0][1] * A[1][0]
    total = 0
    for j in range(n):
        menor = [linha[:j] + linha[j + 1:] for linha in A[1:]]
        total += (-1) ** j * A[0][j] * determinante_laplace(menor)
    return total


# Utilitários de matrizes
def multiplicar(A, B):
    n, m, p = len(A), len(B), len(B[0])
    C = [[0.0] * p for _ in range(n)]
    for i in range(n):
        for k in range(m):
            aik = A[i][k]
            if aik != 0:
                Bk = B[k]
                Ci = C[i]
                for j in range(p):
                    Ci[j] += aik * Bk[j]
    return C


def transposta(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def identidade(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def matriz_aleatoria(n, lo=-9, hi=9):
    return [[random.randint(lo, hi) for _ in range(n)] for _ in range(n)]


def erro_inversa(A, inv):
    """Maior |elemento| de (A * A^-1 - I)."""
    P = multiplicar(A, inv)
    n = len(A)
    return max(abs(P[i][j] - (1 if i == j else 0))
               for i in range(n) for j in range(n))


def aproximadamente_iguais(x, y):
    return abs(x - y) <= TOL * max(1.0, abs(x), abs(y))


# Entrada e saída
def ler_numero(mensagem):
    while True:
        texto = input(mensagem).strip().replace(",", ".")
        try:
            valor = float(texto)
            return int(valor) if valor.is_integer() else valor
        except ValueError:
            print("  Valor inválido. Digite um número.")


def ler_inteiro(mensagem, minimo, maximo):
    while True:
        texto = input(mensagem).strip()
        if texto.isdigit() and minimo <= int(texto) <= maximo:
            return int(texto)
        print(f"  Digite um inteiro entre {minimo} e {maximo}.")


def ler_matriz(n, nome="A"):
    print(f"\nDigite os elementos da matriz {nome} ({n}x{n}):")
    M = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            M[i][j] = ler_numero(f"  {nome}[{i + 1}][{j + 1}] = ")
    return M


def fmt(x):
    if isinstance(x, int):
        return str(x)
    if abs(x) < 1e-12:
        return "0"
    if abs(x - round(x)) < 1e-9 * max(1.0, abs(x)) and abs(x) < 2 ** 53:
        return str(int(round(x)))
    if abs(x) >= 1e9 or abs(x) < 1e-4:
        return f"{x:.6e}"
    return f"{x:.4f}".rstrip("0").rstrip(".")


def mostrar(nome, M):
    largura = max(len(fmt(x)) for linha in M for x in linha)
    print(f"{nome} =")
    for linha in M:
        print("   | " + "  ".join(fmt(x).rjust(largura) for x in linha) + " |")


def titulo(texto):
    print("\n" + "=" * 60)
    print(texto)
    print("=" * 60)


# AÇÕES DO MENU

def acao_ordem_12():
    n = ler_inteiro("\nOrdem da matriz (1 ou 2): ", 1, 2)
    A = ler_matriz(n)
    print()
    mostrar("A", A)
    det = determinante_12(A)
    if n == 1:
        print(f"\ndet(A) = {fmt(det)}")
    else:
        print(f"\ndet(A) = {fmt(A[0][0])}*{fmt(A[1][1])} - "
              f"{fmt(A[0][1])}*{fmt(A[1][0])} = {fmt(det)}")
    inv = inversa_12(A)
    if inv is None:
        print("det(A) = 0 -> matriz singular, não possui inversa.")
    else:
        print()
        mostrar("A^-1", inv)


def acao_sarrus():
    A = ler_matriz(3)
    print()
    mostrar("A", A)
    sarrus_3x3(A)
    print(f"\nConferência (eliminação de Gauss): {fmt(determinante_n(A))}")


def acao_aleatoria():
    n = ler_inteiro("\nOrdem n da matriz aleatória (1 a 200): ", 1, 200)
    A = matriz_aleatoria(n)
    print()
    if n <= 6:
        mostrar("A", A)
    else:
        print(f"(matriz {n}x{n} aleatória, inteiros de -9 a 9; omitida por tamanho)")
    det = determinante_n(A)
    print(f"\ndet(A) = {fmt(det)}")
    inv = inversa_n(A)
    if inv is None:
        print("Matriz singular (det = 0): não possui inversa.")
        return
    if n <= 6:
        print()
        mostrar("A^-1", inv)
    else:
        print("A^-1 calculada (omitida por tamanho).")
    print(f"Verificação: max|A*A^-1 - I| = {erro_inversa(A, inv):.2e}")


def acao_testes():
    random.seed(42)
    titulo("MATRIZES DE TESTE (n = 2, 3, 10, 50, 100)")

    print("\n-- Ordens menores: inteiros aleatórios --")
    for n in (2, 3, 10):
        A = matriz_aleatoria(n)
        det = determinante_n(A)
        inv = inversa_n(A)
        print(f"\nn = {n}: det = {fmt(det)}", end="")
        if n == 2:
            print(f"   (fórmula direta: {fmt(determinante_12(A))})", end="")
        if n == 3:
            print(f"   (Sarrus: {fmt(sarrus_3x3(A, False))})", end="")
        print()
        if n <= 3:
            mostrar("A", A)
        if inv is None:
            print("  singular: sem inversa")
        else:
            print(f"  max|A*A^-1 - I| = {erro_inversa(A, inv):.2e}")

    print("\n-- Ordens maiores: matrizes estruturadas (det conhecido) --")
    casos = [
        (50, "diagonal com todos os elementos = 2", 2, 2 ** 50),
        (100, "matriz identidade", 1, 1),
    ]
    for n, descricao, valor, det_exato in casos:
        A = [[valor if i == j else 0 for j in range(n)] for i in range(n)]
        det = determinante_n(A)
        inv = inversa_n(A)
        erro_rel = abs(det - det_exato) / abs(det_exato)
        print(f"\nn = {n} ({descricao})")
        print(f"  det calculado = {fmt(det)}")
        print(f"  det analítico = {fmt(det_exato)}")
        print(f"  erro relativo = {erro_rel:.2e}")
        print(f"  max|A*A^-1 - I| = {erro_inversa(A, inv):.2e}")


def acao_propriedades():
    n = 5
    titulo("VERIFICAÇÃO DE PROPRIEDADES (n = 5)")

    def verificar(A, B, alpha, linha):
        dA, dB = determinante_n(A), determinante_n(B)
        # 1) Binet
        d_ab = determinante_n(multiplicar(A, B))
        binet = aproximadamente_iguais(d_ab, dA * dB)
        # 2) Transposição
        d_t = determinante_n(transposta(A))
        transp = aproximadamente_iguais(dA, d_t)
        # 3) Multiplicar uma linha por alpha
        C = [l[:] for l in A]
        C[linha] = [alpha * x for x in C[linha]]
        d_c = determinante_n(C)
        escalar = aproximadamente_iguais(d_c, alpha * dA)
        return (dA, dB, d_ab, d_t, d_c, binet, transp, escalar)

    random.seed(7)
    A, B = matriz_aleatoria(n), matriz_aleatoria(n)
    alpha, linha = 3.5, 2
    dA, dB, d_ab, d_t, d_c, binet, transp, escalar = verificar(A, B, alpha, linha)

    mostrar("A", A)
    mostrar("B", B)
    ok = lambda b: "OK" if b else "FALHOU"
    print(f"\ndet(A) = {fmt(dA)}    det(B) = {fmt(dB)}")
    print(f"\n1) Teorema de Binet: det(A*B) = det(A)*det(B)")
    print(f"   det(A*B) = {fmt(d_ab)}   |   det(A)*det(B) = {fmt(dA * dB)}   -> {ok(binet)}")
    print(f"\n2) Transposição: det(A) = det(A^T)")
    print(f"   det(A) = {fmt(dA)}   |   det(A^T) = {fmt(d_t)}   -> {ok(transp)}")
    print(f"\n3) Linha {linha + 1} de A multiplicada por alpha = {alpha}")
    print(f"   det(nova) = {fmt(d_c)}   |   alpha*det(A) = {fmt(alpha * dA)}   -> {ok(escalar)}")

    # Repetição em massa para reforçar a validação
    total = 200
    acertos = [0, 0, 0]
    for _ in range(total):
        X, Y = matriz_aleatoria(n), matriz_aleatoria(n)
        r = verificar(X, Y, random.choice([2, 3.5, -1.5, 0.25]), random.randrange(n))
        for k in range(3):
            acertos[k] += r[5 + k]
    print(f"\nRepetindo com {total} pares aleatórios (tolerância relativa {TOL}):")
    print(f"   Binet: {acertos[0]}/{total}   Transposição: {acertos[1]}/{total}"
          f"   Linha*alpha: {acertos[2]}/{total}")


def _melhor_tempo(f, repeticoes):
    melhor = float("inf")
    for _ in range(repeticoes):
        t0 = time.perf_counter()
        f()
        melhor = min(melhor, time.perf_counter() - t0)
    return melhor


def acao_desempenho():
    titulo("ANÁLISE DE DESEMPENHO")
    random.seed(1)

    print("\n-- Expansão de Laplace (recursiva, O(n!)) x Eliminação de Gauss O(n^3) --")
    print(f"{'n':>4} {'Laplace (s)':>14} {'Gauss (s)':>14} {'Laplace/Gauss':>15}")
    for n in range(2, 10):
        A = [[random.random() for _ in range(n)] for _ in range(n)]
        tl = _melhor_tempo(lambda: determinante_laplace(A), 1)
        tg = _melhor_tempo(lambda: determinante_n(A), 5)
        print(f"{n:>4} {tl:>14.6f} {tg:>14.6f} {tl / tg:>15.1f}")

    print("\n-- Eliminação de Gauss: tempo x ordem --")
    ordens = [10, 20, 40, 80, 160, 240]
    tempos = []
    print(f"{'n':>5} {'tempo (s)':>12} {'t/n^3 (ns)':>12} {'t(n)/t(n ant.)':>16} {'(n/n ant.)^3':>14}")
    for idx, n in enumerate(ordens):
        A = [[random.random() for _ in range(n)] for _ in range(n)]
        t = _melhor_tempo(lambda: determinante_n(A), 3)
        tempos.append(t)
        if idx == 0:
            razao, esperada = "-", "-"
        else:
            razao = f"{t / tempos[-2]:.2f}"
            esperada = f"{(n / ordens[idx - 1]) ** 3:.2f}"
        print(f"{n:>5} {t:>12.6f} {t / n ** 3 * 1e9:>12.1f} {razao:>16} {esperada:>14}")

    expoente = math.log(tempos[-1] / tempos[0]) / math.log(ordens[-1] / ordens[0])
    # expoente usando só as ordens maiores
    expoente2 = math.log(tempos[-1] / tempos[2]) / math.log(ordens[-1] / ordens[2])
    print(f"\nExpoente estimado (inclinação log-log), n=10..240 : {expoente:.2f}")
    print(f"Expoente estimado (inclinação log-log), n=40..240 : {expoente2:.2f}")
    print("Valor teórico esperado: 3 (complexidade O(n^3)).")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(6, 4))
        plt.loglog(ordens, tempos, "o-", label="medido")
        ref = [tempos[2] * (n / ordens[2]) ** 3 for n in ordens]
        plt.loglog(ordens, ref, "--", label="referência O(n³)")
        plt.xlabel("Ordem n")
        plt.ylabel("Tempo (s)")
        plt.title("Determinante por eliminação de Gauss")
        plt.legend()
        plt.grid(True, which="both", alpha=0.3)
        plt.tight_layout()
        plt.savefig("desempenho.png", dpi=150)
        print("Gráfico salvo em desempenho.png")
    except ImportError:
        print("(matplotlib não instalado: gráfico não gerado; use a tabela acima.)")


# MENU
MENU = """
==================== MENU ====================
 1 - Determinante e inversa (ordem 1 ou 2)
 2 - Regra de Sarrus passo a passo (3x3)
 3 - Matriz aleatória de ordem n: det e inversa
 4 - Matrizes de teste (n = 2, 3, 10, 50, 100)
 5 - Verificação de propriedades (n = 5)
 6 - Análise de desempenho e complexidade
 0 - Sair
=============================================="""

ACOES = {
    "1": acao_ordem_12,
    "2": acao_sarrus,
    "3": acao_aleatoria,
    "4": acao_testes,
    "5": acao_propriedades,
    "6": acao_desempenho,
}


def main():
    print("Lista 6 - Determinantes e inversas (implementação manual)")
    while True:
        print(MENU)
        opcao = input("Escolha uma opção: ").strip()
        if opcao == "0":
            print("Encerrando.")
            break
        elif opcao in ACOES:
            ACOES[opcao]()
        else:
            print("Opção inválida.")


if __name__ == "__main__":
    main()
