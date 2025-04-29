import random

import numpy as np
import procedures as proc
import format

def hamming_weight_test1():
    a1 = np.random.randint(2, size=8)
    a2 = np.random.randint(3, size=8)
    a3 = np.random.randint(5, size=8)
    a4 = np.random.randint(3, size=4)
    print("Array: ", a1, ". Hamming weight = ", proc.hamming_weight(a1), ". Is balanced = ", proc.is_balanced(a1))
    print("Array: ", a2, ". Hamming weight = ", proc.hamming_weight(a2), ". Is balanced = ", proc.is_balanced(a2))
    print("Array: ", a3, ". Hamming weight = ", proc.hamming_weight(a3), ". Is balanced = ", proc.is_balanced(a3))
    print("Array: ", a4, ". Hamming weight = ", proc.hamming_weight(a4), ". Is balanced = ", proc.is_balanced(a4))


def hamming_weight_test2():
    a1 = np.random.randint(2, size=16)
    a2 = np.random.randint(2, size=16)
    a3 = np.random.randint(3, size=8)
    a4 = np.random.randint(3, size=8)
    a5 = np.random.randint(2, size=4)
    a6 = np.random.randint(2, size=4)
    print("f : ", a1, "\ng : ", a2, ". Hamming distance = ", proc.hamming_weight(a1, a2))
    print("f : ", a3, "\ng : ", a4, ". Hamming distance = ", proc.hamming_weight(a3, a4))
    print("f : ", a5, "\ng : ", a6, ". Hamming distance = ", proc.hamming_weight(a5, a6))


def format_test():
    f = np.random.randint(2, size=8)
    print("Vector of values :", f)
    pdnf = format.vv_to_pdnf(f)
    print("PDNF :", pdnf)
    print("Backward VV :", format.pdnf_to_vv(pdnf))
    print("ANF :", format.vv_to_anf(f))
    print("ANF to VV:", format.anf_to_vv(format.vv_to_anf(f)))


def tensor_test():
    A = np.random.randint(random.randint(2, 10), size=(random.randint(1, 3), random.randint(1, 3)))
    B = np.random.randint(random.randint(2, 10), size=(random.randint(1, 3), random.randint(1, 3)))
    print("A:")
    for row in A:
        print(row)
    print("B:")
    for row in B:
        print(row)
    C = proc.tensor_product(A, B)
    print("A*B:")
    for row in C:
        print(row)
    C = proc.tensor_product(B, A)
    print("B*A:")
    for row in C:
        print(row)


def mat_mul_test():
    A = np.random.randint(random.randint(2, 10), size=(random.randint(1, 3), random.randint(1, 3)))
    B = np.random.randint(random.randint(2, 10), size=(random.randint(1, 3), random.randint(1, 3)))
    print("A:")
    for row in A:
        print(row)
    print("B:")
    for row in B:
        print(row)
    C = proc.mat_mul(A, B)
    if C is not None:
        print("A*B:")
        for row in C:
            print(row)
    C = proc.mat_mul(B, A)
    if C is not None:
        print("B*A:")
        for row in C:
            print(row)


def wht_test():
    n = 3
    f = np.array([0, 0, 0, 1, 0, 1, 0, 0])
    print("f: ", f)
    w = proc.wht(f)
    print("w : ", w)
    f = np.array([1, 1, 1, 1, 1, 0, 0, 1])
    print("f: ", f)
    w = proc.wht(f)
    print("w : ", w)
    f = np.random.randint(2, size=8)
    print("f: ", f)
    w1 = proc.wht(f)
    print("w1 : ", w1)
    print("f from w1 : ", proc.inv_wht(w1))
    w2 = proc.walsh(f, n)
    print("w2 : ", w2)
    print("f from w2 : ", proc.inv_walsh(w2, n))


def cross_correlation_test():
    f = np.array([0, 1, 1, 0])
    g = np.array([1, 1, 1, 0])
    print("f : ", f)
    print("g : ", g)
    print("Delta(f,g) from def: ", proc.cross_correlation_def(f, g))
    print("Delta(f,f) from def: ", proc.cross_correlation_def(f, f))
    print("Delta(g,g) from def: ", proc.cross_correlation_def(g, g))
    print("Delta(f,g) from theorem: ", proc.cross_correlation(f, g))
    print("Delta(f,f) from theorem: ", proc.cross_correlation(f, f))
    print("Delta(g,g) from theorem: ", proc.cross_correlation(g, g))
    f = np.random.randint(2, size=8)
    g = np.random.randint(2, size=8)
    print("f : ", f)
    print("g : ", g)
    print("Delta(f,g) from def: ", proc.cross_correlation_def(f, g))
    print("Delta(f,f) from def: ", proc.cross_correlation_def(f, f))
    print("Delta(g,g) from def: ", proc.cross_correlation_def(g, g))
    print("Delta(f,g) from theorem: ", proc.cross_correlation(f, g))
    print("Delta(f,f) from theorem: ", proc.cross_correlation(f, f))
    print("Delta(g,g) from theorem: ", proc.cross_correlation(g, g))


def cross_correlation_k_test():
    # Пример полностью некоррелированных функций: x1 и x2
    # f = [0, 1, 0, 1] x1
    # g = [0, 0, 1, 1] x2
    # f = [0, 1, 0, 1, 0, 1, 1, 0] # x1*-x3 + (x1^x2)*x3
    # g = [0, 0, 1, 1, 0, 0, 0, 0] # x2*-x3
    # Порядка k = 1
    # f = [1, 1, 1, 1, 1, 0, 1, 1] # 1 + x1*x3 + x1*x2*x3
    # g = [1, 1, 0, 0, 1, 1, 0, 1] # 1 + x2 + x1*x2*x3
    k = 1  # - Порядок некоррелированности
    f = np.random.randint(2, size=8)
    g = np.random.randint(2, size=8)
    while proc.cross_correlation_k(f, g) != k:
        f = np.random.randint(2, size=8)
        g = np.random.randint(2, size=8)
    print("f : ", f)
    print("g : ", g)
    print("Delta(f,g): ", proc.cross_correlation(f, g), " k : ", proc.cross_correlation_k(f, g))
    print("ANF of f : ", format.vv_to_anf(f))
    print("ANF of g : ", format.vv_to_anf(g))


# Нахождение функции с корр-имунностью равной m
def correlative_immunity_test():
    m = 1
    f = np.random.randint(2, size=8)
    while proc.correlative_immunity_K2(f) != m:
        f = np.random.randint(2, size=8)
    print("f : ", f)
    print("w_f: ", proc.wht(f), " m : ", proc.correlative_immunity_K2(f))
    print("ANF of f : ", format.vv_to_anf(f))


def m_balanced_test():
    f = np.random.randint(2, size=16)
    print(f)
    print("m_balnced m :", proc.m_balanced(f))


# Проверка утверждения 3К - w_f(a) = 0 mod 2^{m+1} если cor(f) = m
def statement_3K():
    m = 0
    # Неравновероятная
    f = np.random.randint(2, size=8)
    t = 0
    while (proc.correlative_immunity_K2(f) != m or proc.is_balanced(f)) and t < 5000:
        f = np.random.randint(2, size=8)
        t += 1
    if t != 5000:
        print("f : ", f)
        print("w_f: ", proc.wht(f), " m : ", proc.correlative_immunity_K2(f), " is balanced: ", proc.is_balanced(f))
        print("ANF of f : ", format.vv_to_anf(f))
        print("Is all w_f(a) = 0 mod 2^{m+2}? :", "YES" if np.sum(proc.wht(f) % pow(2, m + 2)) == 0 else "NO")
    else:
        print("Can't find")
    # Равновероятная
    t = 0
    while (proc.correlative_immunity_K2(f) != m or not proc.is_balanced(f)) and t < 5000:
        f = np.random.randint(2, size=8)
        t += 1
    if t != 5000:
        print("f : ", f)
        print("w_f: ", proc.wht(f), " m : ", proc.correlative_immunity_K2(f), " is balanced: ", proc.is_balanced(f))
        print("ANF of f : ", format.vv_to_anf(f))
        print("Is all w_f(a) = 0 mod 2^{m+2}? :", "YES" if np.sum(proc.wht(f) % pow(2, m + 2)) == 0 else "NO")
    else:
        print("Can't find")


# КР Вариант 1
# f - аффинная функция <B, x> xor 1, B = (1, 1, 0, 0, 1, 1, 0, 0)
# Найти w_f - коэф. Уолша Адамара

def kr_3():
    f1 = np.array([0, 1, 1, 0, 0, 0, 1, 1])
    print(f1, format.vv_to_anf(f1))
    f2 = np.array([1, 0, 0, 1, 0, 1, 1, 0])
    print(f2, format.vv_to_anf(f2))
    print('Cor f1', proc.correlative_immunity_K2(f1))
    print('Cor f2', proc.correlative_immunity_K2(f2))


def BDZ():
    # 514.3
    print('514.3')
    n = 2
    for a in range(1 << (n + 1)):
        f_anf = []
        if a % 2 == 1:
            f_anf.append('1')
        for i in range(1, n + 1):
            if (a >> i) % 2 == 1:
                f_anf.append(f'x_{i}')
        print('f = ', format.anf_to_vv(f_anf, n), 'cor f = ', proc.correlative_immunity_K2(format.anf_to_vv(f_anf, n)))
        # print('f_anf = ', f_anf, 'f = ', format.anf_to_vv(f_anf, n), 'f_anf = ', proc.wht(format.anf_to_vv(f_anf, n)))
    # 515
    # 523
    print('523')
    w_f = np.array([2, -2, 2, 2])
    print('а) w_f =', w_f, ' f = ', proc.inv_wht(w_f), 'f_anf =', format.vv_to_anf(proc.inv_wht(w_f)))
    w_f = np.array([2, 2, -2, -2])
    print('б) w_f =', w_f, ' f = ', proc.inv_wht(w_f), 'f_anf =', format.vv_to_anf(proc.inv_wht(w_f)),
          '- нет такой функции')
    w_f = np.array([-2, 2, 2, 2])
    print('в) w_f =', w_f, ' f = ', proc.inv_wht(w_f), 'f_anf =', format.vv_to_anf(proc.inv_wht(w_f)))


# Перебор всех функций от 3 переменных и вывод их спектров
def perebor():
    n = 3
    for a in range(1 << (1 << n)):
        f = format.a_to_v(a, 1 << n)
        print('w_f = ', proc.wht(f), end='\t')
        print('f = ', f, end='\t')
        print('hamming weight = ', proc.hamming_weight(f), end='\t')
        print('anf:', format.vv_to_anf(f))


def deg_test():
    n = 2
    for a in range(1 << (1 << n)):
        f = format.a_to_v(a, 1 << n)
        print('f = ', f)
        print('anf:', format.vv_to_anf(f))
        print('deg = ', proc.deg(f))


def ann_test():
    # 575
    n = 4
    f = np.zeros(n, dtype=np.dtype(int))
    print("f = ", f, ": ann = ", proc.Ann(f, True))
    f = np.ones(n, dtype=np.dtype(int))
    print("f = ", f, ": ann = ", proc.Ann(f, True))
    f = ['x_1']
    vv = format.anf_to_vv(f, 2)
    print("f = ", vv, ": ann = ", proc.Ann(vv, True))
    f = [['x_1', 'x_2']]
    vv = format.anf_to_vv(f)
    print("f = ", vv, ": ann = ", proc.Ann(vv, True))


def algebraic_immunity_test():
    # f = ['1', ['x_1', 'x_2', 'x_4'], ['x_1', 'x_2']]
    # print(proc.algebraic_immunity_2(format.anf_to_vv(f)))
    # exit()
    f = ['1', ['x_1', 'x_2', 'x_4'], ['x_1', 'x_2']]
    print("f = ", format.anf_to_vv(f), "AI = ", proc.algebraic_immunity_1(format.anf_to_vv(f)))
    f = [['x_1', 'x_2'], ['x_3', 'x_4', 'x_5', 'x_6']]
    not_f = np.array([1 if v == 0 else 0 for v in format.anf_to_vv(f)])
    ai_f = proc.algebraic_immunity_2(format.anf_to_vv(f))
    ai_not_f = proc.algebraic_immunity_2(not_f)
    ai = ai_f if ai_f[0] < ai_not_f[0] else ai_not_f
    print("f = ", format.anf_to_vv(f), "AI = ", ai)


def DDT_and_LAT_test():
    s = [12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2] # s-box of PRESENT
    print("S-box ", s, "\nDDT: ")
    print(proc.DDT(s, 4, 4)) # from GF(2^4) to GF(2^4)
    print()
    print("LAT: ")
    print(proc.LAT(s, 4, 4)) # from GF(2^4) to GF(2^4)
    # Проверка через Sagemath:
    # from sage.crypto.sboxes import PRESENT as s
    # print(s)
    # print(s.difference_distribution_table())
    # print()
    # print(s.linear_approximation_table())
