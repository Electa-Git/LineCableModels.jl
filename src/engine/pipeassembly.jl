function _pipe_surface_terms(wall,position,::Type{T})::Tuple{Complex{T},Complex{T}} where {T}
    outside=wall(Val(:outer))
    return outside,position>1 ? wall(Val(:mutual)) : zero(outside)
end

function _pipe_pair(input,pipe,left,right)
    p=pipe.conductor
    i=last(input.assemblies[left]); j=last(input.assemblies[right])
    centre=input.positions[p]
    zi=(input.positions[i][1]-centre[1],input.positions[i][2]-centre[2])
    zj=(input.positions[j][1]-centre[1],input.positions[j][2]-centre[2])
    radii=(max(input.r_ext[i],input.r_ins_ext[i]),max(input.r_ext[j],input.r_ins_ext[j]))
    return PipeImpedance.Pair(left,right,(zi,zj),radii)
end

function _pipe_interactions(input,methods,rho,s)
    isempty(input.pipes) && return nothing
    return map(input.pipes) do pipe
        i=pipe.conductor
        wall=methods.pipe_impedance(input.r_in[i],input.r_ext[i],rho[i],input.mu_cond[i],s)
        wall.state.proximity===:none && return wall
        rows=map(pipe.children) do child
            length(input.assemblies[child])==1 ||
                throw(ArgumentError("pipe core proximity requires screenless one-conductor coaxial units"))
            row=only(input.assemblies[child])
            iszero(input.r_in[row]) || throw(ArgumentError("pipe core proximity requires solid round cores"))
            row
        end
        positions=[(input.positions[k][1]-input.positions[i][1],
            input.positions[k][2]-input.positions[i][2]) for k in rows]
        cores=PipeImpedance.Cores(pipe.children,positions,input.r_ext[rows],rho[rows],input.mu_cond[rows])
        return PipeImpedance.with_cores(wall,cores)
    end
end

function _pipe_impedance!(destination,input,interactions)
    isempty(input.pipes) && return destination
    groups,_=_assembly_groups(input)
    for (pipe,wall) in zip(input.pipes,interactions)
        p=pipe.conductor
        assembly=findfirst(range->p in range,input.assemblies)
        a=destination[p,p]; transfer=wall(Val(:mutual))
        for (ki,left) in pairs(pipe.children), right in pipe.children[ki:end]
            pair=_pipe_pair(input,pipe,left,right)
            kind=left==right ? Val(:self) : Val(:mutual)
            coefficient=wall(kind,pair)+a-2transfer
            for row in groups[left],column in groups[right]
                destination[row,column]+=coefficient
                left==right || (destination[column,row]+=coefficient)
            end
        end
        for child in pipe.children,row in groups[child],outer in input.assemblies[assembly]
            coefficient=destination[p,outer]-(outer==p ? transfer : zero(transfer))
            destination[row,outer]+=coefficient
            destination[outer,row]+=coefficient
        end
    end
    return destination
end

function _pipe_potential!(destination,input,methods,frequency,temperature,s)
    isempty(input.pipes) && return destination
    groups,_=_assembly_groups(input)
    for pipe in input.pipes
        p=pipe.conductor
        assembly=findfirst(range->p in range,input.assemblies)
        κ=constitutive(methods.insulation_admittance,pipe.material,frequency,temperature)
        iszero(κ) && throw(DomainError(κ,"common pipe cavity has zero admittivity"))
        ε0=one(frequency)*88541878128*(one(frequency)*10)^(-22)
        ε=ε0*pipe.material.eps_r
        cavity=methods.pipe_admittance(input.r_in[p],ε)
        factor=s*ε/κ
        exterior=destination[p,p]
        for (ki,left) in pairs(pipe.children),right in pipe.children[ki:end]
            pair=_pipe_pair(input,pipe,left,right)
            kind=left==right ? Val(:self) : Val(:mutual)
            coefficient=factor*cavity(kind,pair)+exterior
            for row in groups[left],column in groups[right]
                destination[row,column]+=coefficient
                left==right || (destination[column,row]+=coefficient)
            end
        end
        for child in pipe.children,row in groups[child],outer in input.assemblies[assembly]
            coefficient=destination[p,outer]
            destination[row,outer]+=coefficient
            destination[outer,row]+=coefficient
        end
    end
    return destination
end
