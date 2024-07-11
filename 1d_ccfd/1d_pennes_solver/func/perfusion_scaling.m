function val = perfusion_scaling(v_in,v_coma,v_star)

	if (v_in - v_coma) > 0
		val = (v_in - v_coma) / (v_star - v_coma);
	else
		val = 0;
	end

end
