function [partial_x, partial_y] = partial(p)
p_sptr_L = D2Lenphy_specShen(p, 1);
p_x_sptr = LDTL(p_sptr_L);
p_y_sptr = TLLDY(p_sptr_L);
partial_x = D2Lenphy_specShen(p_x_sptr, 2);
partial_y = D2Lenphy_specShen(p_y_sptr, 2);

end