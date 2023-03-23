# Solutions of tasks for "Special Methods of Programming" by D

**Variant №2**.


# Task 1: solve equation
**Task:** find solution of equation:

$$ f(x) = \sin(x) - 1/x = 0 .$$

**Solution:** [here](./task1_solve_equation/src/main.rs).

**Answers:**
```
by dichotomy   : x = 1.114156723022461  , f_evals = 36
by chords      : x = 1.1141571408717854 , f_evals = 24
by newton      : x = 1.114157140871924  , f_evals = 12
by direct iters: x = 1.1141570186671306 , f_evals = 10
```


# Task 1 flex: solve equation at compile time
**Task:** solve "Task 1" at compile time.

**Solution:** [here](./task1_solve_equation_at_compile_time/src/main.rs).

**Answer:**
```
by dichotomy: 1.114156723022461
```


# Task 2: solve polynomial equation
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


# Task 3: find min of Rosenbrock function with fixed precision
**Task:** find min of Rosenbrock function
with fixed precision of $10^{-3}$ in $x$ and $y$.
Initial point is $x=-1.7$, $y=1.7$.

Rosenbrock function:

$$ f(x,y) = (1-x)^2 + 100 (y-x^2)^2 .$$

**Solution:** [here](./task3_find_min_with_fixed_precision/src/main.rs).

**Answers:**
```
by coordinate descent: x = 0.9995093629951177 , y = 0.9990616715937082 , f_evals = 13506
by "fastest"  descent: x = 1.0004988735118054 , y = 1.0009999995954986 , f_evals = 6727068
by downhill   simplex: x = 1.0004704501887318 , y = 1.000914952298626  , f_evals = 194
```


# Task 4: find min of function
**Task:** find min of function:

$$ [1 + (x+y+1)^2 (19-14x+3x^2-14y+6xy+3y^2)] * [30 + (2x-3y)^2 (18-32x+12x^2+48y-36xy+27y^2)] .$$

**Solution:** [here](./task4_find_min/src/main.rs).

**Answer:**
```
by downhill simplex: x = -0.0004520396141742822 , y = -1.0000038171601773
```


# Task 5: Fit data
**Task:** Fit data given in files with some function, find parameters and fit residue.

Functional used as calculate fit residue:

$$ res(f) = \sum_{n=1}^{N} (f(x_n) - y_n)^2 .$$

Algorithm used to minimize fit residue function:
[Pattern search (Hooke-Jeeves method)](https://en.wikipedia.org/wiki/Pattern_search_(optimization)).

Program algorithm:
- Basic version (fit points by given function):
  1. Build function from string (e.g. $f(x)=ax+b$, where $a=1, b=2$).
  2. Setup parameters for pattern search: $\alpha=2.0$, $\beta=1/\alpha$.
  3. Minimize fit residue using pattern search (change parameters).
- Advanced version (automatically build function to fit points by):
  1. Generate random function with given constraints: function complexity, number of parameters, etc.
  2. Fit data by this function.
  3. If residue is less than residue of best function — set this as best.
  4. Repeat forever.

**Solution:** [here](./task5_fit_data/src/main.rs).

## Task 5.1
### Intended fit:
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.1_fit1.png)
### Better fit:
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.1_fit2.png)
### Fun fit:
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.1_fit3.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.1_fit3_zoomed_out.png)

## Task 5.3
### Intended fit:
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.3_fit1.png)
### Better fits:
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.3_fit2.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.3_fit3.png)

## Task 5.4
### Intended fit:
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.4_fit1.png)
### Alternative fits:
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/cfb9e135b18ec9b8ff2ba02e6a87f8f2878a21a0/SMoPD_task5.4_fit2.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/cfb9e135b18ec9b8ff2ba02e6a87f8f2878a21a0/SMoPD_task5.4_fit3.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/cfb9e135b18ec9b8ff2ba02e6a87f8f2878a21a0/SMoPD_task5.4_fit4.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/cfb9e135b18ec9b8ff2ba02e6a87f8f2878a21a0/SMoPD_task5.4_fit5.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/759ff1e26fd01609751deb448dbc5c4bb151c1c8/SMoPD_task5.4_fit6.png)

## Task 5.fun
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit1.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit2.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit3.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit4.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit5.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit6.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit7.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit8.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit9.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit10.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit11.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit12.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit13.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit14.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit15.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit16.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit17.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit18.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit19.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit20.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit21.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit22.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit23.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit24.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit25.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit26.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit27.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit28.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit29.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit30.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit31.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit32.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit33.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit34.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit35.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit36.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit37.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit38.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit39.png)
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/8950b23a09dba13ff678dcceccf518fd47f89ed7/SMoPD_task5.fun_fit40.png)


# Task 6: deconvolution

**Task:** make deconvolution of given data with given instrumental function.

**Solution:** [here](./task6_deconvolution/src/main.rs).

**Answers:**

## Task 6.0:
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/94a9bb994e348f8510ffdab542c1369aa397c263/SMoPD_task6n.0.png)

## Task 6.2:
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/94a9bb994e348f8510ffdab542c1369aa397c263/SMoPD_task6n.2.png)

## Task 6.3:
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/94a9bb994e348f8510ffdab542c1369aa397c263/SMoPD_task6n.3_unfixed.png)

## Task 6.3 (fixed instrumental function):
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/94a9bb994e348f8510ffdab542c1369aa397c263/SMoPD_task6n.3_fixed.png)

## Task 6.4:
![Screenshot](https://raw.githubusercontent.com/dmyTRUEk/images/94a9bb994e348f8510ffdab542c1369aa397c263/SMoPD_task6n.4.png)

