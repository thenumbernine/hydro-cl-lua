__kernel void calcDisplayVar_name(
#if defined(calcDisplayVar_output_tex) 
	__write_only calcDisplayVar_dstImage_t tex,
#elif defined(calcDisplayVar_output_buffer)
	__global real* dest,
#endif
	int displayVar,
	const __global real* buf
) {
	SETBOUNDS(0,0);
	real value = 0;
	int intindex = dim * index;	//side 0
	if (displayVar >= display_wave_0 && displayVar < display_wave_0 + numWaves) {
		const __global real* wave = buf + intindex * numWaves;
		value = wave[displayVar - display_wave_0];
#if 0
	} else if (displayVar >= display_eigen_0 && displayVar < display_eigen_0 + numEigen) {
		const __global real* eigen = buf + intindex * numEigen;
		value = eigen[displayVar - display_eigen_0];
#endif
	} else if (displayVar >= display_deltaUTilde_0 && displayVar < display_deltaUTilde_0 + numWaves) {
		const __global real* deltaUTilde = buf + intindex * numWaves;
		value = deltaUTilde[displayVar - display_deltaUTilde_0];
	} else if (displayVar >= display_rTilde_0 && displayVar < display_rTilde_0 + numWaves) {
		const __global real* rTilde = buf + intindex * numWaves;
		value = rTilde[displayVar - display_rTilde_0];
	} else if (displayVar >= display_flux_0 && displayVar < display_flux_0 + numStates) {
		const __global real* flux = buf + intindex * numStates;
		value = flux[displayVar - display_flux_0];
	} else if (displayVar >= display_deriv_0 && displayVar < display_deriv_0 + numStates) {
		const __global real* deriv = buf + index * numStates;
		value = deriv[displayVar - display_deriv_0];
	} else if (displayVar == display_dt_0) {
		value = buf[index];
	} else if (displayVar == display_orthoError_0) {
		value = buf[intindex];
	} else {
		value = calcDisplayVar_UBuf(displayVar, buf + numStates * index);
	}
#if defined(calcDisplayVar_output_tex) 
	write_imagef(tex, calcDisplayVar_writeImageArgs, (float4)(value, 0., 0., 0.));
#elif defined(calcDisplayVar_output_buffer)
	dest[index] = value;
#endif
}

