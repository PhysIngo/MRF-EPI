function m_out = relx3D(m_in, TE_in, T1_in, T2_in)
% relax vecot by T1 T2
e1te=exp(-TE_in./T1_in);        % relaxation terms
e2te =exp(-TE_in./T2_in);

m_out = zeros(size(m_in));
m_out(1,:) = m_in(1,:).*e2te;
m_out(2,:) = m_in(2,:).*e2te;
m_out(3,:) = m_in(3,:).*e1te + 1 - e1te ;