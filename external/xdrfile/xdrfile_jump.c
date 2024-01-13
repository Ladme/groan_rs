/* Released under MIT License.
 * Copyright (c) 2023-2024 Ladislav Bartos
 *
 * Implementation of functions for jumping in the xdr trajectory files.
 * This file is NOT part of the `xdrfile` library!
 */

#include "xdrfile_jump.h"

/*** PRIVATE FUNCTIONS  ***/

/*! 
 * \brief Get the number of bytes required to jump to the next frame in an xtc file.
 *
 * \param 	xdp				pointer to the XDRFILE
 * \param	target_time		time that is to be reached (in ps)
 * \param	time_precision	smallest distinguishable unit of time
 * 
 * \return	-1 in case of an error
 * 			0 if the target time has been reached or crossed
 * 			number of bytes to reach the start of the next frame
 * 
 * \note	This function was written for the `groan_rs` library.
 */
static int xtc_get_jump_info(XDRFILE *xdp, float target_time, float time_precision) 
{
	// read magic number
	int magic_number = 0;
	if (xdrfile_read_int(&magic_number, 1, xdp) != 1) return -1;

	// check that the magic number is correct
	if (magic_number != 1995) return -1;

	// jump to time information and read it
    if (xdr_jump(xdp, 8) != 0) return -1;
	float time = 0.0;
	if (xdrfile_read_float(&time, 1, xdp) != 1) return -1;

	// check whether the time is higher than or equal to the target time
	if (time >= target_time - time_precision) return 0;

	// get the number of bytes to the next frame
    if (xdr_jump(xdp, 72) != 0) return -1;
	int size = 0;
	if (xdrfile_read_int(&size, 1, xdp) != 1) return -1;
	// add padding
	if (size % 4 != 0) {
		size += 4 - (size % 4);
	}

	return size;
}

/*! 
 * \brief Get the number of bytes required to jump to the next frame in an trr file.
 *
 * \param 	xdp				pointer to the XDRFILE
 * \param	target_time		time that is to be reached (in ps)
 * \param	time_precision	smallest distinguishable unit of time
 * \param   header_size     size of the header in bytes 
 *                          (set by the function if `target_time` is reached)
 * 
 * \return	-1 in case of an error
 * 			0 if the target time has been reached or crossed
 * 			number of bytes to reach the start of the next frame
 * 
 * \note	This function was written for the `groan_rs` library.
 */
static int trr_get_jump_info(XDRFILE *xdp, float target_time, float time_precision, long *header_size) 
{
    // read the entire header
    t_trnheader header = { 0 };
    if (do_trnheader(xdp, 1, &header) != exdrOK) return -1;

    // check time
    if (header.td >= target_time - time_precision) {
        *header_size = 84;
		// if the file is double precision, we add some more bytes
        if (header.bDouble) *header_size += 8;
        return 0;
    }

    // size of the frame should be the sum of the sizes of the individual segments
	// not all sizes are probably currently used; ignoring them could increase the speed of jumping 
	// but it might break the library for some trr files...
    int size = header.ir_size + header.e_size + header.box_size + header.vir_size + header.pres_size + 
    header.top_size + header.sym_size + header.x_size + header.v_size + header.f_size;

    return size;
}

/*** PUBLIC FUNCTIONS  ***/

int xtc_jump_to_start(XDRFILE *xdp, float target_time)
{
	int jump = 0;
	// repeatedly jump over frames until target frame is reached
	while ((jump = xtc_get_jump_info(xdp, target_time, 0.001)) > 0) {
        if (xdr_jump(xdp, jump) != 0) return 1;
	}

	if (jump == -1) return 1;

	// jump back to the start of the frame
    if (xdr_jump(xdp, -16) != 0) return 1;

	return 0;
}

int trr_jump_to_start(XDRFILE *xdp, float target_time)
{
    int jump = 0;
    // set when target is reached and we need to jump back to the start of the frame
    long back = 0;
    // repeatedly jump over frames until target frame is reached
    while ((jump = trr_get_jump_info(xdp, target_time, 0.001, &back)) > 0) {
        if (xdr_jump(xdp, jump) != 0) return 1;
	}

	if (jump == -1) return 1;

	// jump back to the start of the frame
    if (xdr_jump(xdp, -back) != 0) return 1;

	return 0;
}

int xtc_skip_frame(XDRFILE *xdp)
{
	// read magic number
	int magic_number = 0;
	if (xdrfile_read_int(&magic_number, 1, xdp) != 1) {
		// assuming this only fails if we have reached the end of the file
		return 2;
	}
		
	// check that the magic number is correct (check the validity of the frame)
	if (magic_number != 1995) return 1;

	// get the number of bytes to the next frame
	if (xdr_jump(xdp, 84) != 0) return 1;
	int size = 0;
	if (xdrfile_read_int(&size, 1, xdp) != 1) return 1;
	// add padding
	if (size % 4 != 0) {
		size += 4 - (size % 4);
	}

	// this should only fail if we have reached the end of the file 
	// (but does not have to fail)
	if (xdr_jump(xdp, size) != 0) return 2;

	// read magic number (to check whether we are at EOF)
	if (xdrfile_read_int(&magic_number, 1, xdp) != 1) {
		return 2;
	}

	// check validity of the jump
	if (magic_number != 1995) return 1;

	// move back four bytes (to the start of the frame)
	if (xdr_jump(xdp, -4) != 0) return 1;

	return 0;
}

int trr_skip_frame(XDRFILE *xdp)
{
	// read the header of the frame
    t_trnheader header = { 0 };
	// assuming this only fails if we have reached the end of the file
    if (do_trnheader(xdp, 1, &header) != exdrOK) return 2;

    // size of the frame should be the sum of the sizes of the individual segments
    int size = header.ir_size + header.e_size + header.box_size + header.vir_size + header.pres_size + 
    header.top_size + header.sym_size + header.x_size + header.v_size + header.f_size;

    // this should only fail if we have reached the end of the file
	// (but does not have to fail)
	if (xdr_jump(xdp, size) != 0) return 2;

	// try to read the magic number to see whether we are at the end of file
	int magic_number = 0;
	if (xdrfile_read_int(&magic_number, 1, xdp) != 1) return 2;

	// check validity of the jump
	if (magic_number != 1993) return 1;

	// move back four bytes (to the start of the frame)
	if (xdr_jump(xdp, -4) != 0) return 1;

	return 0;
}