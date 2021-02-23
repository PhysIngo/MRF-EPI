function m_out = rx_flip3D(m_in, r_al)
% rotate vectorby an angle arount the x-achsi
% in:   angle in rad
%       phase in rad - currently not working!!! needs to be implented!!!

m_out = zeros(size(m_in));
m_out(1,:) = m_in(1,:);
m_out(2,:) = m_in(2,:).*cos(r_al) + m_in(3,:).*sin(r_al);
m_out(3,:) = - m_in(2,:).*sin(r_al) + m_in(3,:).*cos(r_al);
