using Statistics, StatsBase
using Distributions
using ProgressMeter
using DataFrames
using CSV
using Plots, UnicodePlots

"""
calcu_cor_true(b1, se1, b2, se2; theta = 0.0)

Compute Rb correlation with jackknife SE and p-value.

# Arguments
- `b1` : AbstractVector of betas from dataset 1
- `se1`: AbstractVector of standard errors of betas from dataset 1
- `b2` : AbstractVector of betas from dataset 2
- `se2`: AbstractVector of standard errors of betas from dataset 2
- `theta` : Sample overlap (scalar or vector). Default = 0.

# Returns
NamedTuple with:
- `rb`   : Correlation estimate
- `se_r` : Jackknife SE
- `p`    : p-value
"""
function calcu_cor_true(b1::AbstractVector, se1::AbstractVector,
                        b2::AbstractVector, se2::AbstractVector;
                        theta::Union{Real,AbstractVector} = 0.0)

    n = length(b1)
    @assert n == length(se1) == length(b2) == length(se2) "Input vectors must be same length"

    # Expand scalar theta
    if length(theta) == 1
        theta = fill(theta, n)
    end

    # remove infinities
    mask = .!(isinf.(b1) .| isinf.(b2) .| isinf.(se1) .| isinf.(se2))
    b1, se1, b2, se2, theta = b1[mask], se1[mask], b2[mask], se2[mask], theta[mask]
    n = length(b1)

    # Precompute
    mean_se1sq = mean(se1 .^ 2)
    mean_se2sq = mean(se2 .^ 2)

    var_b1 = var(b1) - mean_se1sq
    var_b2 = var(b2) - mean_se2sq
    var_b1 < 0 && (var_b1 = var(b1))
    var_b2 < 0 && (var_b2 = var(b2))

    cov_b1b2 = cov(b1, b2) - mean(theta) * sqrt(mean_se1sq * mean_se2sq)
    r = cov_b1b2 / sqrt(var_b1 * var_b2)

    # jackknife storage
    r_jack = Vector{Float64}(undef, n)

    # Precompute sums for leave-one-out
    sum_b1, sum_b2 = sum(b1), sum(b2)
    sum_b1sq, sum_b2sq = sum(abs2, b1), sum(abs2, b2)
    sum_b1b2 = sum(b1 .* b2)
    sum_se1sq, sum_se2sq = sum(se1 .^ 2), sum(se2 .^ 2)
    sum_theta = sum(theta)

     for k in 1:n
        n1 = n - 1
        #println(k)
        #flush(stdout)
        mean_b1 = (sum_b1 - b1[k]) / n1
        mean_b2 = (sum_b2 - b2[k]) / n1

        var_b1_jack = ((sum_b1sq - b1[k]^2) / n1 - mean_b1^2) - (sum_se1sq - se1[k]^2) / n1
        var_b1_jack < 0 && (var_b1_jack = (sum_b1sq - b1[k]^2) / n1 - mean_b1^2)

        var_b2_jack = ((sum_b2sq - b2[k]^2) / n1 - mean_b2^2) - (sum_se2sq - se2[k]^2) / n1
        var_b2_jack < 0 && (var_b2_jack = (sum_b2sq - b2[k]^2) / n1 - mean_b2^2)

        cov_b1b2_jack = ((sum_b1b2 - b1[k]*b2[k]) / n1 - mean_b1*mean_b2) -
                        ((sum_theta - theta[k]) / n1 *
                         sqrt(((sum_se1sq - se1[k]^2) / n1) * ((sum_se2sq - se2[k]^2) / n1)))

        r_jack[k] = cov_b1b2_jack / sqrt(var_b1_jack * var_b2_jack)
    end

    r_mean = mean(r_jack)
    se_r = sqrt((n - 1) / n * sum((r_jack .- r_mean).^2))

    stat = (r / se_r)^2
    pval = 1 - cdf(Chisq(1), stat)

    return (rb = r, se_r = se_r, p = pval)
end

df = DataFrame(CSV.File("rb_results/replications_formatted_for_rb_ClumpedVariants_2025-12-17.txt"))

# cis and trans
res=DataFrame(dataset = "NaN", type = "NaN", rb = NaN, rb_se = NaN, p = NaN)
res=res[Not(1), :]

for dataset in unique(df.variable)
    for eQTL_type in ("cis", "trans")
        println("Analysing: " * dataset * " " * eQTL_type * "...")
        df_f = df[(df.variable .== dataset) .& (df.type .== eQTL_type), :]
        df_f = dropmissing(df_f)

        # set theta based on dataset name
        theta_val = occursin("gtex", lowercase(dataset)) ? 1 : 0

        @time begin
            res_temp = calcu_cor_true(df_f.beta, df_f.se, df_f.rep_beta, df_f.rep_se; theta = theta_val)
            res_temp = DataFrame(dataset = dataset, type = eQTL_type, rb = res_temp.rb, rb_se = res_temp.se_r, p = res_temp.p)
        end
           append!(res, res_temp)
           println("Analysing: " * dataset * " " * eQTL_type * "...done!")
    end
end

CSV.write("rb_results/cis_trans_replication_rb_ClumpedVariants_2025-12-17.tsv", res; delim='\t')

# Calculate Rb separately for trans-eQTLs that colocalise with any cell type GWAS and for those that do not

# Read in classified eQTL results
ct = DataFrame(CSV.File("Data/eQtlGwasColocLbfBased/eQTLs_ColocalisingBloodTraitGwas_2025-09-30.txt.gz"))

ct.comb3 = string.(ct.variant, "-", ct.phenotype)
ct_eQTLs_set = Set(ct.comb3)

df2 = df[df.type .== "trans", :]

df2.trans_type = ifelse.(in.(df2.eQTL, Ref(ct_eQTLs_set)), "ct", "intracellular")

res2 = DataFrame(dataset = "NaN", type = "NaN", rb = NaN, rb_se = NaN, p = NaN)
res2 = res2[Not(1), :]

for dataset in unique(df2.variable)
    for eQTL_type in ("ct", "intracellular")
        println("Analysing: " * dataset * " " * eQTL_type * "...")
        df_f = df2[(df2.variable .== dataset) .& (df2.trans_type .== eQTL_type), :]
        df_f = dropmissing(df_f)

        # set theta based on dataset name
        theta_val = occursin("gtex", lowercase(dataset)) ? 1 : 0

        @time begin
            res_temp = calcu_cor_true(df_f.beta, df_f.se, df_f.rep_beta, df_f.rep_se; theta = theta_val)
            res_temp = DataFrame(dataset = dataset, type = eQTL_type, rb = res_temp.rb, rb_se = res_temp.se_r, p = res_temp.p)
        end
           append!(res2, res_temp)
           println("Analysing: " * dataset * " " * eQTL_type * "...done!")
    end
end

CSV.write("rb_results/ccTrans_intracellularTrans_replication_rb_ClumpedVariants_2025-12-17.tsv", res2; delim='\t')

# Count number of variants, genes and eQTLs per comparison

df2.type = df2.trans_type
df2 = select!(df2, Not(:trans_type))

df_count = dropmissing(vcat(df, df2))

transform!(df_count, :eQTL => ByRow(x -> split(x, "-")) => [:variant, :gene])

df_count_grouped = groupby(df_count, [:variable, :type])

res_count_eqtl = combine(df_count_grouped, [:eQTL, :variant, :gene] =>
((a, b, c) -> (a = length(unique(a)),  b = length(unique(b)),  c = length(unique(c)))) =>
[:nr_eQTLs, :nr_variants, :nr_genes])

res_count_eqtl = res_count_eqtl[res_count_eqtl.type .!= "intermediate", :]

CSV.write("rb_results/eQTL_counts_rb_ClumpedVariants_2025-12-17.tsv", res_count_eqtl; delim='\t')
