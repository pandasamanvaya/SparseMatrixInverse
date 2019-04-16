using Test
using IncompleteLU
using LinearAlgebra
using SparseArrays

@testset "Inverse of sparse matrices" begin
    
    #Base case
	A = sprand(2,2, 0.2) + 2I
    @test abs(norm(sparse_mat_inv(A)) - norm(inv(Array(A)))) <= 1e-3
	
    A = sprand(10, 10, 0.2) + 10I
    @test abs(norm(sparse_mat_inv(A)) - norm(inv(Array(A)))) <= 1e-3

    #Checking for dense matrix
    A = sprand(12, 12, 0.8) + 10I
    @test_throws ArgumentError sparse_mat_inv(A)
end