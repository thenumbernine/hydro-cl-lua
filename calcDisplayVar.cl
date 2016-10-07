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
	if (displayVar >= displayFirst_wave && displayVar <= displayLast_wave) {
		const __global real* wave = buf + intindex * numWaves;
		value = wave[displayVar - displayFirst_wave];
	} else if (displayVar >= displayFirst_deltaUTilde && displayVar <= displayLast_deltaUTilde) {
		const __global real* deltaUTilde = buf + intindex * numWaves;
		value = deltaUTilde[displayVar - displayFirst_deltaUTilde];
	} else if (displayVar >= displayFirst_rTilde && displayVar <= displayLast_rTilde) {
		const __global real* rTilde = buf + intindex * numWaves;
		value = rTilde[displayVar - displayFirst_rTilde];
	} else if (displayVar >= displayFirst_flux && displayVar <= displayLast_flux) {
		const __global real* flux = buf + intindex * numStates;
		value = flux[displayVar - displayFirst_flux];
	} else if (displayVar >= displayFirst_deriv && displayVar <= displayLast_deriv) {
		const __global real* deriv = buf + index * numStates;
		value = deriv[displayVar - displayFirst_deriv];
	} else if (displayVar == display_reduce_0) {
		value = buf[index];
	} else if (displayVar == display_error_ortho) {
		value = ((const __global error_t*)buf)[intindex].ortho;
	} else if (displayVar == display_error_flux) {
		value = ((const __global error_t*)buf)[intindex].flux;
	} else if (displayVar >= displayFirst_eigen && displayVar <= displayLast_eigen) {
		value = eigen_calcDisplayVar(displayVar, (const __global eigen_t*)buf + intindex);
	} else {
		value = calcDisplayVar_UBuf(displayVar, buf + numStates * index);
	}
#if defined(calcDisplayVar_output_tex) 
	write_imagef(tex, calcDisplayVar_writeImageArgs, (float4)(value, 0., 0., 0.));
#elif defined(calcDisplayVar_output_buffer)
	dest[index] = value;
#endif
}

