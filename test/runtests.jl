@static if VERSION >= v"0.7-"
    using Test
else
    using Base.Test
end
using MetaRowEchelon

MetaRowEchelon.@linsys (x,y) begin
    x = c
end

MetaRowEchelon.@linsys (x,y) begin
    x = c
    y = d
end

MetaRowEchelon.@linsys (x,y) begin
    2x = c
end

MetaRowEchelon.@linsys (x,y) begin
    2k * dx * x = c
end

MetaRowEchelon.@linsys (x,y) begin
    -x = c
end

MetaRowEchelon.@linsys (x,y) begin
    2x = c + π
end

MetaRowEchelon.@linsys (x,y) begin
    2x = exp(c)
end

MetaRowEchelon.@linsys (x,y) begin
    x + (-1)*x = c
end

MetaRowEchelon.@linsys (x,y) begin
    a*x + b*x = c-π
    d*y + e*x = f+exp(g)
end

linExprs = MetaRowEchelon.@dump_linsys (x,y) begin
    a*x + (b/c)*y = c-sin(π)
    d*y + e*x = f+exp(g)
end
@show linExprs

# # Test_throws
# @test_throws LoadError begin
#     MetaRowEchelon.@dump_linsys (x,y) begin
#     a*x + (b/2)*y + (b/2)*y + 1 = 1
#     d*y + e*x = f
# end
# end

linExprs = MetaRowEchelon.@dump_linsys (x,y) begin
    a*x + (b/2)*y + (b/2)*y = 1
    d*y + e*x = f
end
@show linExprs

# # Test_throws or test_broken
# @test_throws LoadError begin
#     MetaRowEchelon.@dump_linsys (x,y) begin
#     a*x = 1.05
#     d*y + e*x = f
# end
# @show linExprs
# end

# # Test_throws
# @test_throws LoadError begin
#     MetaRowEchelon.@linsys (x,y) begin
#     2*x + 3*x + a*y + (-1)*z = c-π
#     d*y + e*x + π*z = f+exp(g)
# end
# end

linExprs = MetaRowEchelon.@dump_linsys (x,y,z) begin
    2*x + 3*x + a*y + (-1)*z = c-π
    d*y + e*x + π*z = f+exp(g)
end
@show linExprs

linExprs = MetaRowEchelon.@dump_linsys (x,y,z) begin
    2*x + 3*x = g(r)
    d*y + e*x + π*z = f+exp(g)
end
@show linExprs

exprs = @macroexpand MetaRowEchelon.@rref (x,y) begin
    x + y = a
    y = b
end
@show exprs

exprs = @macroexpand MetaRowEchelon.@rref (x,y) begin
    x + 2y = a
    -y = b
end
@show exprs

exprs = @macroexpand MetaRowEchelon.@rref (x,y) begin
    a*x + b*y = c
    e*x + d*y = f
end
@show exprs

#=let a=2, b=3
    exprs = MetaRowEchelon.@rref (x,y) begin
    x + 2y = a
    -y = b
    end
    t=map(eval, exprs)
    @show t
end=#
