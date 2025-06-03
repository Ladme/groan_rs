/* Released under MIT License.
 * Copyright (c) 2023-2025 Ladislav Bartos
 *
 * Implementation of functions for jumping in the xdr trajectory files.
 * This file is NOT part of the `xdrfile` library!
 */

#ifndef _xdrfile_jump_h
#define _xdrfile_jump_h

#include <stdio.h>

#include "xdrfile.h"
#include "xdrfile_trr.h"

/*!
 * \brief Jump to frame of an xtc file corresponding to the specified or higher time ignoring previous frames.
 *
 * \param 	xdp				pointer to the XDRFILE
 * \param	target_time		time that is to be reached (in ps)
 *
 * \return	zero if successful, nonzero in case of an error
 *
 * \note - If the function returns 0, XDRFILE is now set to read the frame corresponding to the provided or closest higher time.
 * \note - This function was written for the `groan_rs` library.
 */
int xtc_jump_to_start(XDRFILE *xdp, float target_time);

/*!
 * \brief Jump to frame of a trr file corresponding to the specified or higher time ignoring previous frames.
 *
 * \param 	xdp				pointer to the XDRFILE
 * \param	target_time		time that is to be reached (in ps)
 *
 * \return	zero if successful, nonzero in case of an error
 *
 * \note - If the function returns 0, XDRFILE is now set to read the frame corresponding to the provided or closest higher time.
 * \note - Jumping through trr files is much slower than jumping through xtc files as we need to read more information from the header.
 *
 * \note - This function was written for the `groan_rs` library.
 */
int trr_jump_to_start(XDRFILE *xdp, float target_time);

/*!
 * \brief Skip one frame of the xtc trajectory.
 *
 * \param 	xdp				pointer to the XDRFILE
 *
 * \return	zero if successful, one in case of an error, two if the end of file has been reached
 *
 * \note - If the function returns 0, XDRFILE is now set to read the next frame after the skipped frame.
 * \note - This function was written for the `groan_rs` library.
 */
int xtc_skip_frame(XDRFILE *xdp);

/*!
 * \brief Skip one frame of the xtc trajectory and set time of the frame that was skipped.
 *
 * \param 	xdp				pointer to the XDRFILE
 * \param   time            pointer to float for storing the simulation time
 *
 * \return	zero if successful, one in case of an error, two if the end of file has been reached
 *
 * \note - If the function returns 0, XDRFILE is now set to read the next frame after the skipped frame.
 * \note - This function was written for the `groan_rs` library.
 */
int xtc_skip_frame_with_time(XDRFILE *xdp, float *time);

/*!
 * \brief Skip one frame of the trr trajectory.
 *
 * \param 	xdp				pointer to the XDRFILE
 *
 * \return	zero if successful, one in case of an error, two if the end of file has been reached
 *
 * \note - If the function returns 0, XDRFILE is now set to read the next frame after the skipped frame.
 * \note - This function was written for the `groan_rs` library.
 */
int trr_skip_frame(XDRFILE *xdp);

/*!
 * \brief Skip one frame of the trr trajectory and set time of the frame that was skipped.
 *
 * \param 	xdp				pointer to the XDRFILE
 * \param   time            pointer to float for storing the simulation time
 *
 * \return	zero if successful, one in case of an error, two if the end of file has been reached
 *
 * \note - If the function returns 0, XDRFILE is now set to read the next frame after the skipped frame.
 * \note - This function was written for the `groan_rs` library.
 */
int trr_skip_frame_with_time(XDRFILE *xdp, float *time);

#endif /* _xdrfile_jump_h */