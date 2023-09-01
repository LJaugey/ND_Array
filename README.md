# N-D Array

This project contains a header file that defines a n-dimensional array (variadic) template. This definition takes advantage of data locality which improves performance over e.g. nested `std::vectors`.
Mathematical expression are also optimzed using lazy evaluation which is implemented with expression templates. Expressions such as `(A+B)*C` is compiled to a `for` loop containing `(A[i]+B[i])*C[i]`. This can also be used when evaluating specific array element since only the specified element is effectively computed (e.g. `(A+B)(i,j,k)` simply compiles to `(A(i,j,k)+B(i,j,k))`).

The main.cpp is a quick overview of the features and the speed_test.cpp is a demonstration of the speed difference between this template and nested `std::valarray`.

**Note:** Must be compiled with c++20

## Usage

### definition
An array with shape `[Dim_1,Dim_2,Dim_3]` is defined with
```
Array<Dim_1,Dim_2,Dim_3> A;
```

**Note:** The copy constructor used to perform a shallow copy to get better performance, however it relied on non-mandatory copy-elision so this behaviour was changed to avoid issues with different compilers. The code has been adapted to remove almost all copies.

### accessors
Elements can be accessed through `Array::operator()` with the right amount of indices:
```
A(4,3,2) = 5;
std::cout<<A(4,3,2)<<std::endl;
```
```
5
```

or with `operator[]` which returns a slice of the original array at the specified index:
```
A[4] = 5;
std::cout<<A[4]<<std::endl;
```
```
5 5 5
5 5 5
5 5 5
```

### operations
Usual arithmetic operations (`+`,`-`,`*`,`/`) and many other standard mathematical functions are implemented in parallel (element-wise).

### Other
The size of each dimension can be accessed with `Array::size(dim)`, where `dim` is the dimension (starts at 0)
