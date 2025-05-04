# Файл перевода представлений функций: Вектор значений (vector of values), СДНФ (Principle Disjunctive Normal Form (PDNF)), АНФ (Algebraic Normal Form (ANF))
import numpy as np


# Значение аргумента функции к виду терма СДНФ (011 -> -x_3*x_2*x_1)
def a_to_x(a, r, mode=0):
    x = []
    t = 1
    while t <= r:
        if a % 2 == 1:
            x.append(f'x_{t}')
        else:
            if mode == 0:
                x.append(f'-x_{t}')
        a >>= 1
        t += 1
    return x


# Значение в бинарный вектор
def a_to_v(a, n=0):
    return np.array(list(map(int, list(bin(a)[2:].zfill(n)))))


# вектор значений в СДНФ
def vv_to_pdnf(vv):
    pdnf = []
    for i in range(len(vv)):
        if vv[i] == 1:
            pdnf.append(a_to_x(i, np.log2(len(vv))))
    return np.array([np.array(row) for row in pdnf])


# СДНФ в вектор значений
def pdnf_to_vv(pdnf):
    if len(pdnf) == 0:
        print("pdnf_to_vv zero function accured")
        return 0
    vv = [0 for _ in range(1 << len(pdnf[0]))]  # размерность функции - длина терма
    for i in range(len(pdnf)):
        a = 0
        t = 1
        for j in pdnf[i]:
            if j[0] != '-':
                a += t
            t <<= 1
        vv[a] = 1
    return np.array(vv)


# Вычисление полинома Жегалкина (АНФ) через Треугольник Паскаля (Ещё можно через БПФ)
def vv_to_anf(vv):
    pascal_row = vv
    anf = []
    for a in range(len(vv)):
        if pascal_row[0] == 1:
            if a == 0:
                anf.append('1')
            else:
                anf.append(a_to_x(a, np.log2(len(vv)), 1))
        new_pascal_row = [0 for _ in range(len(pascal_row) - 1)]
        for i in range(0, len(pascal_row) - 1):
            new_pascal_row[i] = pascal_row[i] ^ pascal_row[i + 1]
        pascal_row = new_pascal_row
    return anf
    # return np.array([np.array(row) for row in anf])


# Подставляет значение a - big-endian в ANF [1, ['x_2'], ['x_3', 'x_2', 'x_1']]
def sub_a(x, a, n=None):
    if x == '1':
        return 1
    if n is None:
        n = get_n_from_anf(x)
    mp = {}
    i = 0
    while i < n:
        mp[f'x_{i + 1}'] = a % 2
        a >>= 1
        i += 1
    res = 0
    for e in x:
        if e == '1':
            res = res ^ 1
        else:
            if type(e) is list: # mode for anfs [[x_1, x_2], [x_3]]
                f = True
                for el in e:
                    if mp[el] != 1:
                        f = False
                        break
                if f:
                    res = res ^ 1
            else:
                res = 1
                if mp[e] != 1: # mode for monoms x_1*x_2*x_3
                    res = 0
                    break
    return res


def get_n_from_anf(anf):
    n = 0
    if len(anf) == 1 and anf[0] == '1':
        return 0
    for e in anf:
        if e == '1':
            continue
        if type(e) is list:
            for el in e:
                if int(el[el.index('_') + 1:]) > n:
                    n = int(el[el.index('_') + 1:])
        else:
            if int(e[e.index('_') + 1:]) > n:
                n = int(e[e.index('_') + 1:])
    return n


# Вычисление вектора значения функции через АНФ
def anf_to_vv(anf, N=None):
    if N is None:
        n = get_n_from_anf(anf)
    else:
        n = N
    if len(anf) == 0:
        return [0] * (1 << n)
    if n == 0:  # its a constant 1 anf
        if N is None:
            return [1]
        else:
            return [1] * (1 << N)
    vv = [0 for _ in range(1 << n)]
    for a in range(1 << n):
        vv[a] = sub_a(anf, a)
    return np.array(vv)
