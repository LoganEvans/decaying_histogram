from math import log as ln


def max_diff(alpha_fast, alpha_slow):
    a, b = 1 - alpha_fast, 1 - alpha_slow
    return ln(ln(b) / ln(a)) / ln(a / b) - 1

def get_count_diff(alpha_fast, alpha_slow):
    n = max_diff(alpha_fast, alpha_slow)
    x = (1.0 - (1.0 - alpha_fast) ** (n + 1)) / alpha_fast
    y = (1.0 - (1.0 - alpha_slow) ** (n + 1)) / alpha_slow
    return y - x

def find_optimal_decay(alpha_slow):
    max_count_diff = 0
    alpha_fast = alpha_slow + 1e-10
    delta = 0.01
    while delta > 1e-100:
        while True:
            c = []
            for d in [x * delta for x in range(4)]:
                c.append(get_count_diff(alpha_fast + d, alpha_slow))
            if c[1] < c[2] and c[2] < c[3]:
                alpha_fast += delta
            else:
                break
        delta *= 0.1
    return alpha_fast


def main():
    alpha0 = 0.005
    alpha1 = 0.001
    n = max_diff(alpha0, alpha1)
    print n
    print get_count_diff(alpha0, alpha1)
    decay = find_optimal_decay(alpha1)
    print decay
    print get_count_diff(decay, alpha1)
    print max_diff(decay, alpha1)

    #for n in range(397, 403):
    #    #n = max_diff(alpha0, alpha1)
    #    print n
    #    x = (1.0 - (1.0 - alpha0) ** (n + 1))
    #    y = (1.0 - (1.0 - alpha1) ** (n + 1))
    #    print x - y

if __name__ == '__main__':
    main()

