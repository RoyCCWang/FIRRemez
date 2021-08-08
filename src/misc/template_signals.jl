function rectfunc(x::T, x_c::T)::T where T
    if abs(x) <= x_c
        return one(T)
    end

    return zero(T)
end


function raisedcosinefunction(x::T, β, 𝑇)::T where T
    return raisedcosinefunction(x,convert(T,β),convert(T,𝑇))
end

"""
    raisedcosinefunc(x::T, β::T, 𝑇::T)::T
Let low := (1-β)/(2𝑇),
    |x| := abs(x),
    high := (1+β)/(2𝑇).

Returns 1 if |x| ≤ low,
Returns 0.5*(1+cos(π*𝑇/β *(|x|-low))) if low ≤ |x| ≤ high.
Returns 0 otherwise.
"""
function raisedcosinefunction(x::T, β::T, 𝑇::T)::T where T
    two = one(T) + one(T)
    low = (one(T)-β)/(two*𝑇)
    high = (one(T)+β)/(two*𝑇)

    if abs(x) <= low
        return one(T)
    elseif low < abs(x) <= high
        return one(T)/two*( one(T) + cos(π*𝑇/β *(abs(x)-low)) )
    end

    return zero(T)
end


# transition band is cosine. This allows satisfies the power complementary condition.
# cos( (yp-y)/(yp-ys) * π/2 ) for a falling edge. yp is pass band, ys is stop band.
function cosinetransitionlowpassfunction(y::T, yp::T, ys::T)::T where T
    abs_y = abs(y)

    if abs_y <= yp
        return one(T)
    elseif yp < abs_y <= ys
        #return cos( (yp-y)/(ys-yp)*π/2 )
        return cos( (yp-abs_y)/(ys-yp)*π/2 )
    end

    return zero(T)
end

function cosinetransitionhighpassfunction(x::T, xp::T, xs::T)::T where T
    abs_x = abs(x)

    if abs_x <= xs
        return zero(T)
    elseif xs < abs_x <= xp
        #return cos( (xp-x)/(xp-xs)*π/2 )
        return cos( (xp-abs_x)/(xp-xs)*π/2 )
    end

    return one(T)
end


# rs is rising edge's stopband.
# rp is rising edge's passband.
# fp is falling edge's passband.
# fs is falling edge's stopband.
function cosinetransitionbandpassfunction(y::T, rs, rp, fp, fs)::T where T
    abs_y = abs(y)

    # rising edge.
    if abs_y <= rs
        return zero(T)
    elseif rs < abs_y <= rp
        #return cos( (rp-y)/(rp-rs)*π/2 )
        return cos( (rp-abs_y)/(rp-rs)*π/2 )
    end

    # falling edge.
    if rp <= abs_y <= fp
        return one(T)
    elseif fp < abs_y <= fs
        #return cos( (fp-y)/(fs-fp)*π/2 )
        return cos( (fp-abs_y)/(fs-fp)*π/2 )
    end

    return zero(T)
end
