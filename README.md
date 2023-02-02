# Solutions of tasks for "Special Methods of Programming" by D

Variant №2.


## Task 1
Task: find solution of equation: `sin(x) - 1/x = 0`.

Solution: [here](./task1_solve_equation/src/main.rs).

Answers:
```
solution_by_dichotomy: 1.114156723022461
solution_by_chords   : 1.1141571408717854
solution_by_newton   : 1.114157140871924
solution_by_direct_iterations: -2.772604525254792
```


## Task 1: flex
Task: solve "Task 1" at compile time.

Solution: [here](./task1_solve_equation_at_compile_time/src/main.rs).

Answer:
```
X = 1.114156723022461
```


## Task 2
Task: find all solutions of equation: `x⁴ + 2x³ + 7x² + 2x + 7 = 0`.

Solution: [here](./task2_solve_polynomial_equation/src/main.rs).

Answers:
```
0.03564470863974995 + 1.0816682537596316 i
0.03564470863974995 - 1.0816682537596316 i
-1.0356447080166609 + 2.2144580135918535 i
-1.0356447080166609 - 2.2144580135918535 i
```


## Task 3
Task: find min of Rosenbrock function with fixed precision.

Solution: [here](./task3_find_min_with_fixed_precision/src/main.rs).

Answers:
```
by coordinate descent: x = 0.9997922598100459 , y = 0.9995361936770234
by fastest    descent: x = 0.9995008705143117 , y = 0.9990000002451329
```

