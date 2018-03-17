/// 786

/******************************************************************************/

#pragma once

/******************************************************************************/

/* Counts the barcodes in the interleaved FASTQ from stdin
 *
 * Parameters:
 *   known_barcodes_path: path for 10x known barcodes
 *   output_prefix:       prefix used to generate prefix.ema-fcnt and
 *                        prefix.ema-ncnt files used by _ema correct_
 *   max_map_size:        maximum in-memory size of the barcode hashmap
 *                        larger values increase memory consumption of both
 *                        count and correct stage, but speed up correct stage
 */
EXTERNC void count(
	const char *known_barcodes_path,
	const char *output_prefix,
	const size_t max_map_size
);
