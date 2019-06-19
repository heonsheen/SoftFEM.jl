abstract type TestAT end

mutable struct TestT <: TestAT
    x::Float64
    v::Float64            
end

function test_func(testat::TestAT)
    (testat.x, testat.v)
end