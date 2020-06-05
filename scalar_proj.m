function [Sdot] = scalar_proj(vec_a,vec_b,magn)

    if size(vec_a,2)==1 && size(vec_b,2)~=1
        vec_a=repmat(vec_a,1,size(vec_b,2));
    elseif size(vec_b,2)==1 && size(vec_a,2)~=1
        vec_b=repmat(vec_b,1,size(vec_a,2));
    end

    Sdot = dot(vec_a,vec_b)./(magn);

end