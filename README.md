# Solutions of tasks for "Special Methods of Programming" by D

**Variant №2**.


## Task 1: solve equation
**Task:** find solution of equation:

$$ \sin(x) - 1/x = 0 .$$

**Solution:** [here](./task1_solve_equation/src/main.rs).

**Answers:**
```
by dichotomy   : 1.114156723022461
by chords      : 1.1141571408717854
by newton      : 1.114157140871924
by direct iters: 1.1141570186671306
```


## Task 1 flex: solve equation at compile time
**Task:** solve "Task 1" at compile time.

**Solution:** [here](./task1_solve_equation_at_compile_time/src/main.rs).

**Answer:**
```
by dichotomy: 1.114156723022461
```


## Task 2: solve polynomial equation
**Task:** find all solutions of equation:

$$ x^4 + 2x^3 + 7x^2 + 2x + 7 = 0 .$$

**Solution:** [here](./task2_solve_polynomial_equation/src/main.rs).

**Answers:**
```
0.03564470863974995 + 1.0816682537596316 i
0.03564470863974995 - 1.0816682537596316 i
-1.0356447080166609 + 2.2144580135918535 i
-1.0356447080166609 - 2.2144580135918535 i
```


## Task 3: find min of Rosenbrock function with fixed precision
**Task:** find min of Rosenbrock function
with fixed precision of $10^{-3}$ in $x$ and $y$.
Initial point is $x=-1.7$, $y=1.7$.

Rosenbrock function:

$$ f(x,y) = (1-x)^2 + 100 (y-x^2)^2 .$$

**Solution:** [here](./task3_find_min_with_fixed_precision/src/main.rs).

**Answers:**
```
by coordinate descent: x = 0.9995093629951177 , y = 0.9990616715937082 , iters = 13506
by fastest    descent: x = 1.0004988735118054 , y = 1.0009999995954986 , iters = 6727068
by downhill   simplex: x = 1.0004704501887318 , y = 1.000914952298626  , iters = 194
```


## Task 4: find min of function
**Task:** find min of function:

$$ [1 + (x+y+1)^2 (19-14x+3x^2-14y+6xy+3y^2)] * [30 + (2x-3y)^2 (18-32x+12x^2+48y-36xy+27y^2)] .$$

**Solution:** [here](./task4_find_min/src/main.rs).

**Answer:**
```
by downhill simplex: x = -0.0004520396141742822 , y = -1.0000038171601773
```

