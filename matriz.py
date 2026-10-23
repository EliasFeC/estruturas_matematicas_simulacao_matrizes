
# funções
def determinante(A):
    """
    Ordem 1: det(A) = a
    Ordem 2: det(A) = a*d - b*c
    """
    n = len(A)
    if n == 1:
        return A[0][0]
    a, b = A[0]
    c, d = A[1]
    return a * d - b * c


def inversa(A):
    """
    Ordem 1: A^-1 = [[1/a]]
    Ordem 2: A^-1 = (1/det) * [[ d, -b],
                               [-c,  a]]
    Retorna None se det = 0
    """
    det = determinante(A)
    if det == 0:
        return None
    if len(A) == 1:
        return [[1 / det]]
    a, b = A[0]
    c, d = A[1]
    return [[ d / det, -b / det],
            [-c / det,  a / det]]


# entrada/saída
def ler_numero(mensagem):
    """Lê um número do teclado"""
    while True:
        texto = input(mensagem).strip().replace(",", ".")
        try:
            valor = float(texto)
            return int(valor) if valor.is_integer() else valor
        except ValueError:
            print("  Valor inválido. Digite um número.")


def ler_ordem():
    """pede a ordem da matriz"""
    while True:
        texto = input("\nOrdem da matriz (1 ou 2): ").strip()
        if texto in ("1", "2"):
            return int(texto)
        print("  Ordem inválida. Digite 1 ou 2.")


def ler_matriz():
    """pede os elementos da matriz"""
    n = ler_ordem()
    print(f"\nDigite os elementos da matriz A ({n}x{n}):")
    M = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            M[i][j] = ler_numero(f"  A[{i + 1}][{j + 1}] = ")
    return M


def fmt(x):
    if isinstance(x, int):
        return str(x)
    texto = f"{x:.4f}".rstrip("0").rstrip(".")
    return "0" if texto == "-0" else texto


def mostrar(nome, M):
    largura = max(len(fmt(x)) for linha in M for x in linha)
    print(f"{nome} =")
    for linha in M:
        print("   | " + "  ".join(fmt(x).rjust(largura) for x in linha) + " |")


# ações
def acao_determinante(A):
    det = determinante(A)
    if len(A) == 1:
        print(f"\ndet(A) = {fmt(A[0][0])}")
    else:
        print(f"\ndet(A) = {fmt(A[0][0])}*{fmt(A[1][1])} - "
              f"{fmt(A[0][1])}*{fmt(A[1][0])} = {fmt(det)}")


def acao_inversa(A):
    inv = inversa(A)
    if inv is None:
        print("\ndet(A) = 0, não possui inversa.")
        return
    print()
    mostrar("A^-1", inv)


# menu
MENU = """
================ MENU ================
 1 - Calcular determinante
 2 - Calcular inversa
 3 - Digitar nova matriz A
 0 - Sair
======================================="""


def main():
    A = ler_matriz()
    print()
    mostrar("A", A)

    while True:
        print(MENU)
        opcao = input("Escolha uma opção: ").strip()

        if opcao == "0":
            print("Encerrando.")
            break
        elif opcao == "1":
            acao_determinante(A)
        elif opcao == "2":
            acao_inversa(A)
        elif opcao == "3":
            A = ler_matriz()
            print()
            mostrar("A", A)
        else:
            print("Opção inválida.")


if __name__ == "__main__":
    main()