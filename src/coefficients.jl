export Coefficient

#TODO: Keep closer track of the structure of the coefficient to enable
# simple implementation of algebraic simplification routines.
# In particular, TODO: implement Sequential Comparative Simplification
struct Coefficient{T1,T2}
    num::T1
    den::T2
    function Coefficient(num::T1, den::T2) where {T1<:Integer, T2<:Integer}
        @assert den != 0
        if num == 0
            new{T1, Int}(zero(T1), one(T2))
        elseif num == den
            new{T1, T2}(one(T1), one(T2))
        else
            new{T1, T2}(num, den)
        end
    end
    Coefficient(num::T1, den::T2) where {T1, T2} = new{T1, T2}(num, den)
end
Coefficient(num::T) where {T} = Coefficient(num, 1)



# Product
_coeffProd(a::Coefficient, b::Coefficient) = Coefficient(:($(a.num) * $(b.num)), :($(a.den) * $(b.den)))
function _coeffProd(a::Coefficient{T1,<:Integer}, b::Coefficient{<:Integer,T2}) where {T1,T2}
    if a.den == b.num
        Coefficient(a.num, b.den)
    elseif a.den == 1
        Coefficient(:($(a.num) * $(b.num)), b.den)
    elseif b.num == 1
        Coefficient(a.num, :($(a.den) * $(b.den)))
    else
        Coefficient(:($(a.num) * $(b.num)), :($(a.den) * $(b.den)))
    end
end
_coeffProd(a::Coefficient{<:Integer,T1}, b::Coefficient{<:Integer,T2}) where {T1,T2} = Coefficient(a.num * b.num, :($(a.den) * $(b.den)))
_coeffProd(a::Coefficient{T1, <:Integer}, b::Coefficient{T2, <:Integer}) where {T1,T2} = Coefficient(:($(a.num) * $(b.num)), a.den * b.den)

# Non-simplifying rational number arithmetic. TODO: Pass through a rational number to avoid overflow.
coeffProd(a::Coefficient{<:Integer, <:Integer}, b::Coefficient{<:Integer, <:Integer}) = Coefficient(a.num * b.num, a.den * b.den)
coeffProd(a::Coefficient{<:Integer, <:Integer}, b::Coefficient) = (a.num == a.den) ? b : _coeffProd(a, b)
coeffProd(a::Coefficient, b::Coefficient{<:Integer,<:Integer}) = coeffProd(b, a)

coeffProd(a, b) = _coeffProd(a, b)

# Quotient
coeffQuot(a, b) = coeffProd(a, Coefficient(b.den,b.num))


# Fallthrough 
function _coeffSum(a::Coefficient, b::Coefficient)
    if a.den == b.den
        if a.num == b.num
            Coefficient(:(2*$(a.num)), a.den)
        else    
            Coefficient(:($(a.num) + $(b.num)), a.den)
        end
    else
        Coefficient(:($(a.num)*$(b.den) + $(b.num)*$(a.den)), :($(a.den)*$(b.den)))
    end
end
function _coeffSum(a::Coefficient{<:Integer,T}, b::Coefficient{<:Integer, T}) where {T}
    if a.den == b.den
        Coefficient(a.num + b.num, a.den)
    else
        Coefficient(:($(a.num)*$(b.den) + $(b.num)*$(a.den)), :($(a.den)*$(b.den)))
    end
end

# Non-simplifying rational arithmetic
function coeffSum(a::Coefficient{<:Integer,<:Integer}, b::Coefficient{<:Integer, <:Integer})
    if a.den == b.den
        Coefficient(a.num + b.num, a.den)
    else
        Coefficient((a.num)*(b.den) + (b.num)*(a.den), (a.den)*(b.den))
    end
end
coeffSum(a::Coefficient{<:Integer,<:Integer}, b::Coefficient) = (a.num == 0) ? b : _coeffSum(a, b)
coeffSum(a::Coefficient, b::Coefficient{<:Integer,<:Integer}) = coeffSum(b, a)

coeffSum(a, b) = _coeffSum(a, b)



# Fallthrough 
function _coeffDiff(a::Coefficient, b::Coefficient)
    if a.den == b.den
        if a.num == b.num
            Coefficient(0)
        else    
            Coefficient(:($(a.num) - $(b.num)), a.den)
        end
    else
        Coefficient(:($(a.num)*$(b.den) - $(b.num)*$(a.den)), :($(a.den)*$(b.den)))
    end
end
function _coeffDiff(a::Coefficient{T1, <:Integer}, b::Coefficient{T2, <:Integer}) where {T1,T2}
    if a.den == b.den
        if a.num == b.num
            Coefficient(0)
        else    
            Coefficient(:($(a.num) - $(b.num)), a.den)
        end
    else
        if a.den == 1 && b.den == 1
            Coefficient(:($(a.num) - $(b.num)), 1)
        elseif a.den == 1
            Coefficient(:($(a.num)*$(b.den) - $(b.num)), b.den)
        elseif b.den == 1
            Coefficient(:($(a.num) - $(b.num)*$(a.den)), a.den)
        else
            Coefficient(:($(a.num)*$(b.den) - $(b.num)*$(a.den)), (a.den)*(b.den))
        end
    end
end
function _coeffDiff(a::Coefficient{<:Integer,T}, b::Coefficient{<:Integer, T}) where {T}
    if a.den == b.den
        if a.num == b.num
            Coefficient(0)
        else
            Coefficient(a.num - b.num, a.den)
        end
    else
        if a.num == 1 && b.num == 1
            Coefficient(:($(b.den) - $(a.den)), :($(a.den)*$(b.den)))
        elseif a.num == 1
            Coefficient(:($(b.den) - $(b.num)*$(a.den)), :($(a.den)*$(b.den)))
        elseif b.num == 1
            Coefficient(:($(a.num)*$(b.den) - $(a.den)), :($(a.den)*$(b.den)))
        else
            Coefficient(:($(a.num)*$(b.den) - $(b.num)*$(a.den)), :($(a.den)*$(b.den)))
        end
    end
end

# Non-simplifying rational arithmetic
function coeffDiff(a::Coefficient{<:Integer,<:Integer}, b::Coefficient{<:Integer, <:Integer})
    if a.den == b.den
        Coefficient(a.num - b.num, a.den)
    else
        Coefficient((a.num)*(b.den) - (b.num)*(a.den), (a.den)*(b.den))
    end
end

coeffDiff(a::Coefficient{<:Integer,<:Integer}, b::Coefficient) = (a.num == 0) ? coeffProd(Coefficient(-1), b) : _coeffDiff(a, b)
coeffDiff(a::Coefficient, b::Coefficient{<:Integer,<:Integer}) = (b.num == 0) ? a : _coeffDiff(a, b)
coeffDiff(a, b) = _coeffDiff(a, b)


#TODO: Logic is duplicated between simplification and 
# arithmetic operations. Probably defer all simplification
# to the end?

# Output/simplification routines, primarily for pretty-printing.
_extract_factors(a::Symbol) = [1, Dict(a => 1)]
_extract_factors(a::T) where {T <: Integer} = [a, Dict{Union{Expr, Symbol}, Int}()]
function _extract_factors(a::Expr)
    local numFactor = 1
    local factors = Dict{Union{Expr, Symbol}, Int}()
    if (a.head == :call) && (first(a.args) == :(*)) && (length(a.args) == 3)
        # Monomial-type term
        arg1 = a.args[2] # Extract factors from multiplier
        if arg1 isa Symbol
            if haskey(factors, arg1)
                factors[arg1] += 1
            else
                factors[arg1] = 1
            end
        elseif arg1 isa Int
            numFactor *= arg1
        else
            arg1_num, arg1_factors = _extract_factors(arg1)
            for (arg, num) in arg1_factors
                if haskey(factors, arg)
                    factors[arg] += 1
                else
                    factors[arg] = 1
                end
            end
        end
        
        arg2 = a.args[3] # Extract factor from multiplicand
        if arg2 isa Symbol
            if haskey(factors, arg2)
                factors[arg2] += 1
            else
                factors[arg2] = 1
            end
        elseif arg1 isa Int
            numFactor *= arg2
        else
            arg2_num, arg2_factors = _extract_factors(arg2)
            for (arg, num) in arg2_factors
                if haskey(factors, arg)
                    factors[arg] += 1
                else
                    factors[arg] = 1
                end
            end
        end
    # elseif a.head == :call && (first(a.args) == :(+) || first(a.args) == :(-)) && (length(a.args) == 3)
    #     # Handle simple sums or differences
    #     op = first(a.args)
    #     arg1 = a.args[2]
    #     arg2 = a.args[3]
    #     if ((arg1 isa Expr && arg1.head == :(*)) && (arg2 isa Expr && arg2.head == :(*)))
    #         const1, factors1 = _extract_factors(arg1)
    #         const2, factors2 = _extract_factors(arg2)
    #         for term in keys(factors1)
    #             if term in keys(factors2)
    #                 nrepeats = min(factors1[term], factors2[term])
    #                 if haskey(factors, term)
    #                     factors[term] += nrepeats
    #                 else
    #                     factors[term] = nrepeats
    #                 end
    #                 factors1[term] -= nrepeats
    #                 factors2[term] -= nrepeats
    #             end
    #         end
    #         remaining1=:($(const1))
    #         for (factor, ntimes) in factors1
    #             if ntimes == 1
    #                 remaining1 = :($(remaining1)*$(factor))
    #             else
    #                 remaining1 = :($(remaining1)*$(factor)^$(ntimes))
    #             end
    #         end
    #         remaining2=:($(const2))
    #         for (factor, ntimes) in factors2
    #             if ntimes == 1
    #                 remaining2 = :($(remaining2)*$(factor))
    #             else
    #                 remaining1 = :($(remaining2)*$(factor)^$(ntimes))
    #             end
    #         end
    #         remainingTerm = :($(remaining1)$(op)$(remaining2))
    #         if haskey(factors, remainingTerm)
    #             factors[remainingTerm] += 1
    #         else
    #             factors[remainingTerm] = 1
    #         end
    #     else
    #         if haskey(factors, a)
    #             factors[a]+=1
    #         else
    #             factors[a]=1
    #         end
    #     end
    else# Not a product or simple sum of terms, so cannot factor (yet)
        if haskey(factors, a)
            factors[a] += 1
        else
            factors[a] = 1
        end
    end
    numFactor, factors
end
_simplify(a::Union{Symbol,<:Integer}) = a
#TODO: This code (much like the arithmetic code above) begs for some abstraction.
# Let's not ask so much of the compiler to simplify these, terms, particularly 
# since this arithmetic is "exact".
function _simplify(a::Expr)
    if (a.head == :call) && (first(a.args) == :(*)) && (length(a.args) == 3)
        # Product
        arg1 = _simplify(a.args[2])
        arg2 = _simplify(a.args[3])
        if arg1 == 1
            return arg2
        elseif arg2 == 1
            return arg1
        elseif (arg1 == 0 || arg2 == 0)
            return 0
        else
            return :($(arg1)*$(arg2))
        end
    elseif (a.head == :call) && (first(a.args) == :(-)) && (length(a.args) == 3)
        # Difference
        arg1 = _simplify(a.args[2])
        arg2 = _simplify(a.args[3])
        if arg1 == 0
            return _simplify(:((-1)*$(arg2)))
        elseif arg2 == 0
            return arg1
        else
            return :($(arg1)-$(arg2))
        end
    elseif (a.head == :call) && (first(a.args) == :(+)) && (length(a.args) == 3)
        # Sum
        arg1 = _simplify(a.args[2])
        arg2 = _simplify(a.args[3])
        if arg1 == 0
            return arg2
        elseif arg2 == 0
            return arg1
        else
            return :($(arg1)+$(arg2))
        end
    elseif a.head == :call && (first(a.args) == :(/)) && (length(a.args) == 3)
        arg1 = a.args[2]
        num1, arg1_factors = _extract_factors(arg1)
        arg2 = a.args[3]
        num2, arg2_factors = _extract_factors(arg2)
        for a1 in keys(arg1_factors) # Merge factors according to power rules
            for (a2, n2) in arg2_factors
                if a1 == a2
                    if arg1_factors[a1] == n2
                        delete!(arg1_factors, a1)
                        delete!(arg2_factors, a2)
                    elseif arg1_factors[a1] > n2
                        arg1_factors[a1] -= n2
                        delete!(arg2_factors, a1)
                    elseif arg1_factors[a1] < n2
                        arg2_factors[a1] -= arg1_factors[a1]
                        delete!(arg1_factors, a1)
                    end
                end
            end
        end
        arg1_out = num1
        for (term,numTimes) in arg1_factors
            if arg1_out == 1
                if numTimes == 1
                    arg1_out = term
                else
                    arg1_out = :($(term)^$(numTimes))
                end
            else
                if numTimes == 1
                    arg1_out = :($(arg1_out)*$(term))
                else
                    arg1_out = :($(arg1_out)*$(term)^$(numTimes))
                end
            end
        end
        arg2_out = num2
        for (term,numTimes) in arg2_factors
            if arg2_out == 1
                if numTimes == 1
                    arg2_out = term
                else
                    arg2_out = :($(term)^$(numTimes))
                end
            else
                if numTimes == 1
                    arg2_out = :($(arg2_out)*$(term))
                else
                    arg2_out = :($(arg2_out)*$(term)^$(numTimes))
                end
            end
        end
        return :($(arg1_out)/$(arg2_out))
    else
        # @warn("Cannot simplify with head $(a.head) and args $(a.args).")
        return a
    end
end
function _simplify(a::Coefficient) 
    aden_simplify = _simplify(a.den)
    anum_simplify = _simplify(a.num)
    if aden_simplify == 1
        anum_simplify
    elseif anum_simplify == 0
        return 0
    elseif anum_simplify == aden_simplify
        return 1
    else
        return _simplify(:($(anum_simplify) / $(aden_simplify)))
    end
end
function simplify(a::Coefficient{<:Integer, <:Integer})
    if a.den == 1
        a.num
    elseif a.den == -1
        -a.num
    else
        Rational(a.num, a.den)
    end
end
function simplify(a::Coefficient{T, <:Integer}) where {T} 
    if (a.den == 1)
        _simplify(a.num)
    elseif a.den == -1
        _simplify(:($(-1)*$(a.num)))
    else
        _simplify(a)
    end
end
simplify(a::Coefficient{<:Integer, T}) where {T} = (a.num == 0) ? 0 : _simplify(a)
simplify(a::Coefficient) = _simplify(a)