import math
from itertools import combinations
from typing import List, Set, Dict, Any

def get_prime_factorization(n: int) -> List[int]:
    """
    Returns a list of the prime factors of a positive integer.
    """
    factors = []
    while n % 2 == 0:
        factors.append(2)
        n //= 2
    p = 3
    while p * p <= n:
        while n % p == 0:
            factors.append(p)
            n //= p
        p += 2
    if n > 1:
        factors.append(n)
    return factors

def calculate_lcm_from_primes(primes: Set[int]) -> int:
    """Calculates the LCM of a set of distinct primes, which is their product."""
    if not primes:
        return 1
    return math.prod(primes)

def solve_for_partition(m_x: int, m_y: int, n: int, s_x: Set[int], s_y: Set[int]) -> List[Dict[str, Any]]:
    """
    Finds all solutions for x + y = n, where x is a multiple of m_x and y is a multiple of m_y.
    """
    solutions = []
    for j in range(1, (n // m_y) + 1):
        remainder = n - (m_y * j)
        if remainder <= 0:
            break
        if remainder % m_x == 0:
            x = remainder
            y = m_y * j
            solution = {
                "x": x,
                "y": y,
                "S_x (x的质因数)": sorted(list(s_x)),
                "m_x (x的最小公倍数)": m_x,
                "S_y (y的质因数)": sorted(list(s_y)),
                "m_y (y的最小公倍数)": m_y
            }
            solutions.append(solution)
    return solutions

def find_all_solutions_at_least_3(N: int):
    """
    Main function to find all solutions for a given N where x and y are divisible by at least 3 primes from S.

    Parameters
    ----------
    N : int
        The target integer such that x + y = N.
    """
    S = {2, 3, 5, 7, 11, 13, 17}
    all_solutions = []
    
    size_of_s_x = 3
    for s_x_tuple in combinations(S, size_of_s_x):
        s_x_3_primes = set(s_x_tuple)
        s_y_4_primes = S - s_x_3_primes
        
        m_3 = calculate_lcm_from_primes(s_x_3_primes)
        m_4 = calculate_lcm_from_primes(s_y_4_primes)
        
        solutions1 = solve_for_partition(m_3, m_4, N, s_x_3_primes, s_y_4_primes)
        all_solutions.extend(solutions1)
        
        solutions2 = solve_for_partition(m_4, m_3, N, s_y_4_primes, s_x_3_primes)
        all_solutions.extend(solutions2)

    # 去重: x 和 y 可互换, 只保留一组 (x, y) / (y, x)
    unique_solutions: Dict[tuple, Dict[str, Any]] = {}
    for sol in all_solutions:
        # 使用无序键 (较小值, 较大值) 来判重
        key = tuple(sorted((sol["x"], sol["y"])) )
        if key not in unique_solutions:
            unique_solutions[key] = sol

    return list(unique_solutions.values())

# --- 主要逻辑从这里开始 ---
if __name__ == "__main__":
    print("条件: x和y需分别被集合S中至少3个质数整除。")
    # 获取用户输入的 N
    while True:
        try:
            N_input = int(input("请输入一个正整数 N: "))
            if N_input <= 0:
                raise ValueError
            break
        except ValueError:
            print("无效输入，请输入一个大于零的正整数。")

    print("\n步骤1: 正在计算所有可能的解集...")
    solutions_list = find_all_solutions_at_least_3(N_input)
    print(f"计算完成！共找到 {len(solutions_list)} 组解。")
    
    print("\n步骤2: 按 |x-y| 排序, 选出初始前20名。")
    sorted_by_diff = sorted(solutions_list, key=lambda sol: abs(sol['x'] - sol['y']))
    top_20_by_diff = sorted_by_diff[:20]
    
    print("\n步骤3: 对初始前20名进行二次排序，并保留最大质因数最小的前10名...")
    print("二次排序规则: x 或 y 的最大质因数, 数值越小排名越靠前。")

    # --- 新增逻辑: 计算二次排序键 ---
    for sol in top_20_by_diff:
        factors_x = get_prime_factorization(sol['x'])
        factors_y = get_prime_factorization(sol['y'])
        sol['max_prime_factor'] = max(max(factors_x), max(factors_y))
        sol['fx_str'] = " * ".join(map(str, factors_x))
        sol['fy_str'] = " * ".join(map(str, factors_y))
    # ------------------------------------

    # --- 二次排序 ---
    reranked_solutions = sorted(top_20_by_diff, key=lambda sol: sol['max_prime_factor'])

    # 保留前 10 名
    final_solutions = reranked_solutions[:10]

    print("\n步骤4: 输出最终结果（前10名）。")
    print("-" * 80)
    
    print(f"最终排名 (基于最大质因数): \n")
    for i, sol in enumerate(final_solutions):
        max_pf = sol['max_prime_factor']
        
        print(f"最终排名 #{i+1}: (最大质因数 = {max_pf})")
        print(f"  x = {sol['x']:,}")
        print(f"  y = {sol['y']:,}")
        print(f"  -> x 的S集内质因数 S_x = {sol['S_x (x的质因数)']} (个数: {len(sol['S_x (x的质因数)'])})")
        print(f"  -> y 的S集内质因数 S_y = {sol['S_y (y的质因数)']} (个数: {len(sol['S_y (y的质因数)'])})")
        print(f"  -> x 的质因数完全分解: {sol['fx_str']}")
        print(f"  -> y 的质因数完全分解: {sol['fy_str']}")
        print()
            
    print("-" * 80)