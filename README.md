# Solutions of tasks for "Special Methods of Programming" by D

**Variant №2**.


## Task 1: solve equation
**Task**: find solution of equation: $\sin(x) - 1/x = 0$.

**Solution**: [here](./task1_solve_equation/src/main.rs).

**Answers**:
```
by dichotomy: 1.114156723022461
by chords   : 1.1141571408717854
by newton   : 1.114157140871924
by direct_iterations: -2.772604525254792
```


## Task 1 flex: solve equation at compile time
**Task**: solve "Task 1" at compile time.

**Solution**: [here](./task1_solve_equation_at_compile_time/src/main.rs).

**Answer**:
```
by dichotomy: 1.114156723022461
```


## Task 2: solve polynomial equation
**Task**: find all solutions of equation: $x^4 + 2x^3 + 7x^2 + 2x + 7 = 0$.

**Solution**: [here](./task2_solve_polynomial_equation/src/main.rs).

**Answers**:
```
0.03564470863974995 + 1.0816682537596316 i
0.03564470863974995 - 1.0816682537596316 i
-1.0356447080166609 + 2.2144580135918535 i
-1.0356447080166609 - 2.2144580135918535 i
```


## Task 3: find min of Rosenbrock function with fixed precision
**Task**: find min of Rosenbrock function
with fixed precision of $10^{-3}$ in $x$ and $y$.

Rosenbrock function:

$$ f(x,y) = (1-x)^2 + 100 (y-x^2)^2 $$

**Solution**: [here](./task3_find_min_with_fixed_precision/src/main.rs).

**Answers**:
```
by coordinate descent: x = 0.9997922598100459 , y = 0.9995361936770234
by fastest    descent: x = 0.9995008705143117 , y = 0.9990000002451329
by downhill   simplex: x = 0.9995644823869441 , y = 0.9991330575392396
```

