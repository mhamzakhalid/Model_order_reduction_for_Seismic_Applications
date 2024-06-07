function [u_ortho] = Gram_Schmidt(Vk, u, Xh)

if ~isempty(Vk)
    
    alpha = Vk'*(Xh*u);
    
    u_ortho    = u - Vk*alpha;
    
else
    
    u_ortho    = u;
    
end
    u_ortho = u_ortho/sqrt(abs(dot(u_ortho,(Xh*u_ortho))));
    
end