abstract type TestAT end

mutable struct TestT <: TestAT
    x::Float64
    v::Float64            
end

function test_func(testat::TestAT)
    undef
end

function test_func(testat::TestT)
    (testat.x, testat.v)
end

function test_func2(testat::TestAT)
    test_func(testat)
end