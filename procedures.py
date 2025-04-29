import math
import operator

import numpy as np
import galois
import pandas as pd
from tqdm import tqdm

import format
import gf2


# m-равновероятность, сильная равновероятность
def m_balanced(f, m=None):
    n = int(np.log2(len(f)))
    anf_f = format.vv_to_anf(f)
    if m is None:
        m = pow(2, n - 1)
    for i in range(1, pow(2, n - 1)):  # i - количество гамм (m-граммы)
        freq = [0 for _ in range(pow(2, i))]
        mp = {}
        for a in range(0, pow(2, n + i - 1)):
            for x in range(1, n + i):
                mp[f'x_{x}'] = a % 2
                a >>= 1
            m_gramm = 0
            for j in range(0, i):
                result = 0
                for element in anf_f:
                    if element == '1':
                        result = result ^ 1
                        continue
                    f = True
                    for e in element:
                        if mp[f'x_{int(e[2:]) + j}'] != 1:
                            f = False
                            break
                    if f:
                        result = result ^ 1
                m_gramm = m_gramm + (result << j)  # собираем i-грамму
            freq[m_gramm] += 1
        for j in freq:
            if j != pow(2, n - 1):
                return f"{i - 1}-равновероятная"
    return "Сильно равновероятная "


# Степень функции
def deg(anf_f):
    if isinstance(anf_f, np.ndarray):  # Если дали вектор значений
        anf_f = format.vv_to_anf(anf_f)
    r = 0
    for e in anf_f:
        if e != 1:
            r = max(len(e), r)
    return r

# Вычисляет нелинейность булевой функции 
def nonlinearity(f):
    n = int(np.log2(len(f)))
    return pow(2, n - 1) - max(wht(f)) // 2


# Множество аннигиляторов(аннуляторов) функции f (вектор значений)
# @bool_anf - флаг формата уравнений выходного множества
def Ann(f, bool_anf=False):
    ann = []
    wt = hamming_weight(f)
    pos = [0] * (len(f) - wt)
    j = 0
    for i in range(0, len(f)):
        if f[i] == 0:
            pos[j] = i
            j += 1
    for a in range(pow(2, len(f) - wt)):
        t = a
        g = [0] * len(f)
        for p in pos:
            g[p] = t % 2
            t = t >> 1
        if bool_anf:
            ann.append(format.vv_to_anf(g))
        else:
            ann.append(g)
    if len(ann) == 0:
        if bool_anf:
            ann.append([])
        else:
            ann.append(np.ones(n, dtype=np.dtype(int)))
    if bool_anf:
        ann.sort(key=lambda x: (deg(x), len(x)))
    return ann


# Порядок алгебраической иммуности (Вариант через полный перебор)
def algebraic_immunity_1(f):
    not_f = [1 if v == 0 else 0 for v in f]
    ann1 = Ann(f, True)
    ann2 = Ann(not_f, True)
    r = float('inf')
    if len(ann1) > 1:
        r = deg(ann1[1])  # 0 - тривиальный аннулятор
    if len(ann2) > 1:
        r = min(r, deg(ann2[1]))  # 0 - тривиальный аннулятор
    print(ann1[1], ann2[1])
    return r


# Порядок алгебраической иммуности (Вариант через решение слу в gf(2))
def algebraic_immunity_2(f):
    n = int(np.log2(len(f)))
    monoms = {0: '1'}
    a = 1
    wt = 1
    t = 1
    while a <= pow(2, n):
        monom = []
        a_vv = format.a_to_v(a, n)
        for i in range(n - 1, -1, -1):
            if a_vv[i] == 1:
                monom.append(f'x_{n - i}')
        monoms[t] = monom
        t += 1
        if a == pow(2, n) - 1:
            break
        a = next_ham(a)
        if a >= pow(2, n):
            wt += 1
            a = pow(2, wt) - 1
    supp = [0 for _ in range(hamming_weight(f))]
    t = 0
    for i in range(pow(2, n)):
        if f[i] == 1:
            supp[t] = i
            t += 1
    GF2_matrix = None
    front = 0
    solution = []
    r = 0
    for d in range(0, pow(2, n)):  # d = deg of monoms
        new_monoms = [[0 for _ in range(math.comb(n, d))] for _ in range(hamming_weight(f))]
        for j in range(math.comb(n, d)):
            for i in range(hamming_weight(f)):  # supp(f)
                new_monoms[i][j] = format.sub_a(monoms[front + j], supp[i],
                                                n)  # значение мономов (по i) для a из supp(f)
        if GF2_matrix is None:
            GF2_matrix = np.array(new_monoms)
        else:
            GF2_matrix = np.concatenate((GF2_matrix, np.array(new_monoms)), axis=1)
        front += math.comb(n, d)
        solution_gen = gf2.solve_gf2(GF2_matrix, b=[0 for _ in range(hamming_weight(f))])
        solution = gf2.peek(solution_gen)
        if solution is not None:
            if not np.array(solution).any():
                solution = gf2.peek(solution_gen)
            if solution is not None:
                r = d
                break
    result_ann = []
    for i in range(len(solution)):
        if solution[i] == 1:
            result_ann.append(monoms[i])
    return r, result_ann


# DDT (Differential Distribution Table)
def DDT(f, n, m):
    ddt = np.zeros((pow(2, n), pow(2, m)), dtype=int)
    for i in range(pow(2, n)):
        for j in range(pow(2, m)):
            ddt[j ^ i][f[j] ^ f[i]] += 1
    df = pd.DataFrame(ddt)
    # df.to_excel('ddt.xlsx', index=False, header=False)
    return df


# LAT (Linear Approximation Table)
def LAT(f, n, m):
    lat = np.zeros((pow(2, n), pow(2, m)), dtype=int)
    for i in range(pow(2, n)):
        for j in range(pow(2, m)):
            for x in range(pow(2, n)):
                lat[i][j] = lat[i][j] + (scalar_product(i, x) == scalar_product(j, f[x]))
            lat[i][j] -= pow(2, n-1)
    df = pd.DataFrame(lat)
    # df.to_excel('ddt.xlsx', index=False, header=False)
    return df


# Функция вычисления взаимной корреляции между функциями f и g - по определению
def cross_correlation_def(f, g):
    if len(f) != len(g):
        print("Cross_corelation incorrect size of vectors!")
        return None
    res = []
    for e in range(len(f)):
        ans = 0
        for b in range(len(f)):
            ans += pow(-1, f[b] ^ g[b ^ e])
        res.append(ans)
    return np.array(res)


# Функция вычисления взаимной корреляции между функциями f и g - через Теорему о взаимной корреляции
def cross_correlation(f, g):
    wf = wht(f)
    wg = wht(g)
    d = np.frompyfunc(lambda x, y: x * y, 2, 1)(wf, wg)
    # Создание матрицы Адамара - Сильвестра
    H_1 = np.array([[1, 1], [1, -1]])
    H = np.array([[1, 1], [1, -1]])
    for i in range(int(np.log2(len(f))) - 1):
        H = tensor_product(H, H_1)
    # Делим на 2^n - чтобы получить H_n^{-1}
    inv_H = H / len(f)
    v = mat_mul(d, inv_H).astype(np.int32)
    return v


# Подсчёт порядка некоррелированности (k)
def cross_correlation_k(f, g):
    d = cross_correlation(f, g)
    if d[0] != 0:
        return -1
    for k in range(1, int(np.log2(len(f)))):
        a = pow(2, k) - 1
        while a < len(f):
            if d[a] != 0:
                return k - 1
            a = next_ham(a)
    return int(np.log2(len(f)))  # Абсолютно некоррелированные


# Bill Gosper's algorithm to generate next number with the same hamming weight
def next_ham(a):
    c = a & -a
    r = a + c
    b = (((r ^ a) >> 2) // c) | r
    return b


# Порядок корреляционной иммунности функции f через второй критерий(К2)
def correlative_immunity_K2(f):
    wf = wht(f)
    for m in range(1, int(np.log2(len(f))) + 1):
        a = pow(2, m) - 1
        while a < len(f):
            if wf[a] != 0:
                return m - 1
            a = next_ham(a)
    return int(np.log2(len(f)))


# Вес Хемминга (для одной функции) и Расстояние Хемминга (между двумя функцями)
# Есть другой алгоритм - // Reingold and Nievergelt counting hamming weight in O(1)
def hamming_weight(f, g=None):
    if g is None:
        return np.sum(np.frompyfunc(lambda x: 1 if x != 0 else 0, 1, 1)(f))
    else:
        if len(f) != len(g):
            print("Hamming_weight incorrect size of vectors!")
            exit()
        return np.sum(np.frompyfunc(lambda x, y: 0 if x == y else 1, 2, 1)(f, g))


def is_balanced(f):
    return hamming_weight(f) == (len(f) >> 1)


# Скалярное произведение значений a,b или векторов a,b
def scalar_product(a, b):
    res = 0
    if a is not list and not isinstance(a, np.ndarray):
        while a > 0 or b > 0:
            res ^= (a % 2) * (b % 2)  # Значения - двумерные вектора V_m(2)
            a >>= 1
            b >>= 1
    else:
        for i in range(min(len(a), len(b))):
            res = (res + a[i] * b[i]) % n  # Значения - n-мерные вектора V_m(n)
    return res


# Тензорное произведение матриц A и B (A * B)
def tensor_product(A, B):
    q = len(A) * len(B)
    r = len(A[0]) * len(B[0])
    C = [[0 for _ in range(r)] for _ in range(q)]
    for i in range(len(A)):
        for j in range(len(A[0])):
            for k in range(len(B)):
                for s in range(len(B[0])):
                    C[i * len(B) + k][j * len(B[0]) + s] = A[i][j] * B[k][s]
    return np.array(C)


# Произведение матриц A (r x q) и B (q x k)
def mat_mul(A, B):
    if A.ndim == 1:  # Если A - вектор
        A = A.reshape((1, len(A)))
    if B.ndim == 1:  # Если B - вектор
        B = B.reshape((1, len(B)))
    if len(A[0]) != len(B):
        print("mat_mul sizes of matrix are incorrect : ", len(A), 'x', len(A[0]), " and ", len(B), 'x', len(B[0]))
        return None
    C = [[0 for _ in range(len(B[0]))] for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(A[0])):
                C[i][j] = C[i][j] + A[i][k] * B[k][j]
    C = np.array(C)
    if C.shape[0] == 1:  # Если результат - вектор
        C = C.reshape((C.shape[1],))
    return C


# Преобразование Уолша Адамара через матрицу (свои операции)
# Выход - коэффициенты Уолша Адамара 2-го рода
def wht(f):
    # Создание матрицы Адамара - Сильвестра
    H_1 = np.array([[1, 1], [1, -1]])
    H = np.array([[1, 1], [1, -1]])
    for i in range(int(np.log2(len(f))) - 1):
        H = tensor_product(H, H_1)
    v = np.array([-1 if f[i] == 1 else 1 for i in range(len(f))])
    w = mat_mul(v, H)
    return w


# Обратное преобразование
def inv_wht(w):
    # Создание матрицы Адамара - Сильвестра
    H_1 = np.array([[1, 1], [1, -1]])
    H = np.array([[1, 1], [1, -1]])
    for i in range(int(np.log2(len(w))) - 1):
        H = tensor_product(H, H_1)
    # Делим на 2^n - чтобы получить H_n^{-1}
    inv_H = H / len(w)
    v = mat_mul(w, inv_H)
    f = np.array([1 if v[i] == -1 else 0 for i in range(len(v))])
    return f


# По определению
def walsh(f, n):
    res = []
    for a in range(1 << n):
        ans = 0
        for b in range(1 << n):
            ans += pow(-1, f[b] ^ scalar_product(a, b))
        res.append(ans)
    return np.array(res)


# По определению
def inv_walsh(w, n):
    f = []
    for b in range(1 << n):
        sm = 0
        for a in range(1 << n):
            sm += w[a] * pow(-1, scalar_product(a, b))
        if sm > 0:
            f.append(0)
        else:
            f.append(1)
    return f
