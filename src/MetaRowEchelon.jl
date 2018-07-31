module MetaRowEchelon

include("coefficients.jl")

#TODO: Separate into explicit phases with non-exported macros for
# testing purposes:
# 1) Parsing/extracting coefficients of linear equations
# 2) Simplification (?)
# 3) Symbolic rref. Emit information on solvability?

# 
"""
    handle_lhs!(coefficients, vars, lhsExpr, eqNumber)

Given an expression `lhsExpr`, separate the coefficients on the linear terms
`vars` from the coefficients, which are placed into the row `eqNumber` at the
appropriate column.

Note: this routine is v0.6-only. Will require re-implementation for v0.7+.
"""
function handle_lhs! end

function handle_lhs!(coefficients, vars, lhsExpr::Symbol, eqNumber)
    varInd = findfirst(vars, lhsExpr)
    for i in 1:length(vars)
        if i == varInd
            coefficients[eqNumber, i] = Coefficient(1)
        else
            coefficients[eqNumber, i] = Coefficient(0)
        end
    end
end
function handle_lhs_product(vars, exprArgs)
    coeffs = Dict{Symbol,Coefficient}(var => Coefficient(0) for var in vars)
    # Separate a product into a defined variable and everything else.
    remainingArgs = exprArgs
    currVar = :(nothingInternal)
    for arg in remainingArgs
        if arg in vars
            currVar = arg
        end
    end
    if currVar == :(nothingInternal)
        error("None of the variables $(vars) found in term $(exprArgs). Perhaps there's an issue?")
    end
    for arg in remainingArgs
        if arg in vars
            continue
        else
            if coeffs[currVar] == Coefficient(0)
                coeffs[currVar] = Coefficient(arg)
            else
                coeffs[currVar] = coeffProd(coeffs[currVar], Coefficient(arg))
            end
        end
    end
    coeffs
end
function handle_lhs!(coefficients, vars, lhsExpr::Expr, eqNumber)
    coeffs = Dict{Symbol, Coefficient}(var => Coefficient(0) for var in vars)
    if lhsExpr.head == :call && first(lhsExpr.args) == :(*)
        # Product
        coeffs = handle_lhs_product(vars, lhsExpr.args[2:end])
    elseif lhsExpr.head == :call && first(lhsExpr.args) == :(-) && (length(lhsExpr.args) == 2)
        # Negation
        coeffs[last(lhsExpr.args)] = Coefficient(-1)
    elseif lhsExpr.head == :call && first(lhsExpr.args) == :(+) && (length(lhsExpr.args) > 2)
        # Summation
        for term in lhsExpr.args[2:end]
            # Iterate over terms in sum
            if term isa Symbol
                # Term in sum is a "bare" symbol
                if term ∈ vars
                    coeffs[term] = Coefficient(1)
                else
                    dump(lhsExpr)
                    error("Bare, non-variable found on LHS.") 
                end
            elseif term isa Int
                error("Term $(term) on LHS has no variable.")
            elseif term.head == :call && first(term.args) == :(*)
                # Product term
                tcoeff = handle_lhs_product(vars, term.args[2:end])
                for (var, coeff) in tcoeff
                    if coeff != Coefficient(0)
                        coeffs[var] = coeffSum(coeffs[var], coeff)
                    end
                end
            else
                dump(lhsExpr)
                error("LHS of this form not currently supported.")
            end
        end
    else
        dump(lhsExpr)
        error("LHS of this form not currently supported.")
    end

    for (var, coeff) in coeffs
        varInd = findfirst(vars, var)
        coefficients[eqNumber, varInd] = coeff
    end
end

"""
    handle_rhs!(rhsVector, rhsExpr, eqNumber)

Parses the right-hand side `rhsExpr` into a `Coefficient`, doing relatively
little checking, and places the resulting `Coefficient` into `rhsVector`
at position `eqNumber`.
"""
function handle_rhs! end

function handle_rhs!(rhsVector, rhsExpr::Symbol, eqNumber)
    rhsVector[eqNumber] = Coefficient(rhsExpr)
end
function handle_rhs!(rhsVector, rhsExpr::Expr, eqNumber)
    if rhsExpr.head == :block
        for arg in rhsExpr.args
            if arg isa Symbol || arg isa Int
                # What does the appearance of multiple args mean?
                rhsVector[eqNumber] = Coefficient(arg)
            elseif arg.head == :line
                continue # Line number node
            elseif arg.head == :call
                handle_rhs!(rhsVector, arg, eqNumber)
            else
                dump(arg)
                error("RHS of this form (arg) not currently supported.")
            end
        end
    elseif rhsExpr.head == :call
        rhsVector[eqNumber] = Coefficient(rhsExpr)
    else
        dump(rhsExpr)
        error("RHS of this form not currently supported.")
    end
end

"""
    grok_equality!(coefficients, rhsVector, vars, exprs, eqNumber)

Parse an equality already separated into the vector `exprs` as a
left-hand side which is linear with respect to `vars` and a right-
hand side into `rhsVector`. It is assumed that the corresponding
equation number `eqNumber` is known.
"""
function grok_equality!(coefficients, rhsVector, vars, exprs, eqNumber)
    @assert length(exprs) == 2
    lhsExpr = first(exprs)
    handle_lhs!(coefficients, vars, lhsExpr, eqNumber)

    rhsExpr = last(exprs)
    handle_rhs!(rhsVector, rhsExpr, eqNumber)
end

"""
    parse_linsys!(coefficients, rhsVector, vars, codeLines)

Fill the matrix `coefficients` and right-hand size vector `rhsVector`
with the linear factors and right-hand sides, respectively, of the
equations (linear in `vars` given in `codeLines` (as an AST).
"""
function parse_linsys!(coefficients, rhsVector, vars, codeLines)
    eqNumber = 0
    for codeLine in codeLines
        if codeLine.head == :(=)
            eqNumber += 1
            grok_equality!(coefficients, rhsVector, vars, codeLine.args, eqNumber)
        elseif codeLine.head == :line
            continue
        else
            dump(codeLine)
            error("Line unhandled: ", codeLine)
        end
    end

    coefficients, rhsVector, eqNumber
end

"""
    dump_linear_eqs(coefficients, rhsVector, vars, eqNumber)

Given a coefficient matrix `coefficients`, right-hand side vector `rhsVector`,
variables `vars`, and number of equations `eqNumber`, pretty-print the linear
equations represented therein.
"""
function dump_linear_eqs(coefficients, rhsVector, vars, eqNumber)
    exprs = Expr[]
    numVars = length(vars)
    for eq in 1:eqNumber
        currCoeff = simplify(coefficients[eq, 1])
        currRhs = simplify(rhsVector[eq])
        if numVars == 1
            push!(exprs, :($(currCoeff) * $(vars[1]) = $(currRhs)))
        else
            currExpr = :($(currCoeff) * $(vars[1]))

            for j in 2:numVars
                var = vars[j]
                currCoeff = simplify(coefficients[eq, j])
                currExpr = :($(currExpr) + $(currCoeff)*$(var))
            end
            push!(exprs, :($(currExpr) = $(currRhs)))
        end
    end
    exprs
end

"""
    prepare_containers(vars, blockLines)

Accept a raw list of variables `vars` and collection of block lines `blockLines`
and form matrix of coefficients and right-hand side vector of appropriate size.
"""
function prepare_containers(vars, blockLines)
    numVars = length(vars)
    numEq = length(blockLines) # This is an overestimate, since other, non-equation nodes may be present.
    coefficientMatrix = Matrix{Coefficient}(numEq, numVars)
    rhsVector = Vector{Coefficient}(numEq)
    coefficientMatrix, rhsVector
end

"""
    @linsys(vars::Expr, exprs::Expr)

Return the coefficient matrix, right-hand side vector, and number of equations
for the linear system `exprs` in variables `vars`.
"""
macro linsys(vars::Expr, exprs::Expr)
    if exprs.head == :block
        vars = vars.args
        blockLines = exprs.args
        coefficientMatrix, rhsVector = prepare_containers(vars, blockLines)
        coefficients, rhsVector, eqNumber = parse_linsys!(coefficientMatrix, rhsVector, vars, blockLines)
        :($(coefficients), $(rhsVector), $(eqNumber))
    else
        error("linsys currently only works on blocks.",
        "Please file an issue with your use-case.")
    end
end

"""
    @dump_linsys(vars::Expr, exprs::Expr)

Pretty-print the linear system given in `exprs` with respect to variables `vars`.
"""
macro dump_linsys(vars::Expr, exprs::Expr)
    if exprs.head == :block
        vars = vars.args
        blockLines = exprs.args
        coefficientMatrix, rhsVector = prepare_containers(vars, blockLines)
        coefficients, rhsVector, eqNumber = parse_linsys!(coefficientMatrix, rhsVector, vars, blockLines)
        pretty_exprs = dump_linear_eqs(coefficients, rhsVector, vars, eqNumber)
        :($(pretty_exprs))
    else
        error("dump_linsys currently only works on blocks.",
        "Please file an issue with your use-case.")
    end
end

#TODO: Generalize fwd_elim! and bkwd_elim!, since they accomplish nearly the same goal.
function _fwd_elim!(coefficientMatrix, rhsVector, currCol, numEqs)
    if currCol < 1
        return nothing
    end
    for row in 1:(currCol-1)
        rhsVector[row] = coeffDiff(
            rhsVector[row],
            coeffProd(rhsVector[currCol], coefficientMatrix[row, currCol])
            )
    end
    _fwd_elim!(coefficientMatrix, rhsVector, currCol - 1, numEqs)
    return nothing
end
function _bkwd_elim!(coefficientMatrix, rhsVector, currCol, numEqs)
    if currCol == numEqs
        rhsVector[currCol] = coeffQuot(rhsVector[currCol], coefficientMatrix[currCol, currCol])
        coefficientMatrix[currCol, currCol] = Coefficient(1)
        return nothing
    end
    local pivotFound = false
    local pivotRow = 0
    local pivot = Coefficient(0)
    for currRow in currCol:numEqs
        if pivotFound
            beneathPivot = coefficientMatrix[currRow, currCol]
            for col in currCol:numEqs
                thisProd = coeffProd(beneathPivot, coefficientMatrix[pivotRow, col])
                coefficientMatrix[currRow, col] = coeffDiff(
                    coefficientMatrix[currRow, col],
                    thisProd
                    )
            end
            rhsVector[currRow] = coeffDiff(
                    rhsVector[currRow],
                    coeffProd(beneathPivot, rhsVector[pivotRow])
                )
            continue
        end
        pivot = coefficientMatrix[currRow, currCol]
        pivotRow = currRow
        if pivot != Coefficient(0)
            # Pivot found for current column currCol in row currRow
            pivotFound = true
            for col in (currCol+1):numEqs
                coefficientMatrix[currRow, col] = coeffQuot(
                    coefficientMatrix[currRow, col],
                    pivot
                )
            end
            rhsVector[currRow] = coeffQuot(
                rhsVector[currRow],
                pivot
            )
            coefficientMatrix[currRow, currCol] = Coefficient(1)
        end
    end
    if pivotFound == false
        error("Pivot not found in column $(currCol).")
    end
    
    _bkwd_elim!(coefficientMatrix, rhsVector, currCol + 1, numEqs)
end

function rref!(coefficientMatrix, rhsVector, numEqs)
    numRows, numCols = size(coefficientMatrix)
    @assert numEqs == numCols "System is definitely not well-posed, which is not currently supported."
    _bkwd_elim!(coefficientMatrix, rhsVector, 1, numEqs)

    _fwd_elim!(coefficientMatrix, rhsVector, numEqs, numEqs)
end

"""
    @rref(vars::Expr, exprs::Expr)

Solve the linear system given in `exprs` for the variables `vars` using Gaussian elimination.

Note: Assumes that the system is (symbolically) well-posed, in that each variable
has at least one nonzero coefficient.
"""
macro rref(vars::Expr, exprs::Expr)
    if exprs.head == :block
        vars = vars.args
        blockLines = exprs.args
        coefficientMatrix, rhsVector = prepare_containers(vars, blockLines)
        coefficients, rhsVector, eqNumber = parse_linsys!(coefficientMatrix, rhsVector, vars, blockLines)
        rref!(coefficients, rhsVector, eqNumber)
        

        Expr(:block,
        [ :($(var) = $(simplify(rhsVector[i]))) for (i, var) in enumerate(vars)]
        )
    else
        error("rref currently only works on blocks. Please file an issue with your use-case.")
    end
end

end # module
