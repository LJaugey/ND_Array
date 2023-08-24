# N-D Array

This project contains a header file that defines a n-dimensional array (variadic) template. This definition takes advantage of data locality which improves performance over e.g. nested `std::vectors`.
Mathematical expression are also optimzed using lazy evaluation which is implemented with expression templates. Expressions such as `(A+B)*C` is compiled to a `for` loop containing `(A[i]+B[i])*C[i]`. This can also be used when evaluating specific array element since only the specified element is effectively computed (e.g. `(A+B)(i,j,k)` simply compiles to `(A(i,j,k)+B(i,j,k))`).

The main.cpp is a quick overview of the features and the speed_test.cpp is a demonstration of the speed difference between this template and nested `std::valarray`.

**Note:** Must be compiled with c++17

## Usage

### definition
An array with shape `[Dim_1,Dim_2,Dim_3]` is defined with
```
Array<Dim_1,Dim_2,Dim_3> A;
```

**Warning:** To get better performance, the copy constructor only performs a shallow copy. Deep copies must be explicit with `Array::copy()` or with `Array::operator=(Array& other)`

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
Usual arithmetic operations (`+`,`-`,`*`,`/`) are implemented in parallel and operations `min()`,`max()`,`abs()` are also available
**NOTE:** min/max do not work yet with expression templates.

### Other
The size of each dimension can be accessed with `Array::size(dim)`, where `dim` is the dimension (starts at 0)
