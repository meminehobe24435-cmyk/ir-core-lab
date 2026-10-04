/*
 * frame.h -- BT.656 style video timing and the binary serial frame protocol
 * ------------------------------------------------------------------------
 * TIMING
 * ------
 * A classic BT.656 (625 line / 50 Hz) stream is 864 pixels per line of which
 * 720 carry active video, and 625 lines per frame of which 576 are active,
 * with two clock cycles per pixel (13.5 MHz pixel clock on a 27 MHz clock):
 *
 *      clocks / line  = (720 + 144) * 2 = 1728
 *      clocks / frame = 1728 * 625       = 1 080 000   (40 ms at 27 MHz)
 *
 * A 128x128 IR core is represented with the same structure but a much smaller
 * raster (128 active + 32 blanking pixels per line, 128 active + 16 blanking
 * lines, one clock per pixel):
 *
 *      clocks / line  = 160
 *      clocks / frame = 160 * 144 = 23 040
 *
 * The generator is a plain state machine clocked one pixel clock at a time;
 * it reports data enable, horizontal and vertical sync and the position in
 * the raster, and it counts active pixels and lines so that the structure can
 * be verified instead of assumed.
 *
 * Simplifications, stated honestly: there are no SAV/EAV timing reference
 * codes, no F/V/H bit encoding and no 4:2:2 chroma multiplexing.  Only the
 * *raster structure* of BT.656 is reproduced.
 *
 * SERIAL FRAME PROTOCOL
 * ---------------------
 *      | 0xA5 | 0x5A | type | len_lo | len_hi | payload ... | crc_lo | crc_hi |
 *      |<-- 5 bytes header -->|<-- len bytes -->|<-- 2 bytes -->|
 *
 * The CRC is CRC-16/CCITT-FALSE: polynomial 0x1021, init 0xFFFF, no bit
 * reflection, no final XOR, computed over type + length + payload (i.e. the
 * header magic itself is not protected, exactly like most framing schemes).
 * Its standard check vector is "123456789" -> 0x29B1, asserted in the tests.
 */

#ifndef FRAME_H
#define FRAME_H

#include <stddef.h>
#include <stdint.h>

#include "ir_lab.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/* Video timing                                                        */
/* ------------------------------------------------------------------ */

typedef struct {
    int h_active;          /* active pixels per line            */
    int h_blank;           /* blanking pixels per line          */
    int v_active;          /* active lines per frame            */
    int v_blank;           /* blanking lines per frame          */
    int clocks_per_pixel;  /* 1, or 2 for BT.656                */
} frame_timing_cfg_t;

typedef struct {
    long clk;              /* clock index inside the frame          */
    int  x;                /* pixel index inside the line           */
    int  y;                /* line index inside the frame           */
    int  de;               /* data enable: this clock is active video */
    int  active_pixel;     /* first clock of an active pixel        */
    int  hsync;            /* horizontal sync asserted              */
    int  vsync;            /* vertical sync asserted                */
} frame_signal_t;

typedef struct {
    const frame_timing_cfg_t *cfg;
    long total_clocks;
    long line_clocks;
    long frame_clocks_done;
    /* counters for the frame currently being emitted */
    long active_pixels;
    long lines_done;
    long active_lines_done;
    /* counters of the frame that just completed */
    long last_active_pixels;
    long last_lines;
    long last_active_lines;
    int  frames_done;
} frame_timing_t;

typedef struct {
    long total_clocks;        /* clocks the raster should take            */
    long clocks_emitted;      /* clocks actually walked                   */
    long active_pixel_clocks; /* DE asserted clocks (pixels * clocks_per_pixel) */
    long blanking_clocks;     /* total - active                           */
    long lines;
    long active_lines;
    long pixels_emitted;      /* source pixels copied into the stream     */
    long hsync_asserted;
    long vsync_asserted;
    int  dims_match;
    int  counts_ok;           /* 1 when every structural check passed     */
} frame_report_t;

#define FRAME_ERR_DIM (-20)
#define FRAME_ERR_CAP (-21)

void frame_preset_ir128(frame_timing_cfg_t *c);
void frame_preset_bt656_625(frame_timing_cfg_t *c);

long frame_line_clocks(const frame_timing_cfg_t *c);
long frame_total_clocks(const frame_timing_cfg_t *c);

void frame_timing_init(frame_timing_t *st, const frame_timing_cfg_t *c);
/* Advances one pixel clock.  Returns 1 when the frame just completed. */
int  frame_timing_step(frame_timing_t *st, frame_signal_t *sig);

/*
 * Walks one complete frame, copying active pixels out of `img` (img_w x
 * img_h, row major) into `stream` (one byte per clock, 0x00 during blanking).
 * `stream` may be NULL or too small: the raster is still walked and the
 * counters are still produced, only the copy is skipped/truncated.
 */
int frame_emit(const frame_timing_cfg_t *c, const uint8_t *img,
               int img_w, int img_h, uint8_t *stream, size_t cap,
               frame_report_t *rep);

/* ------------------------------------------------------------------ */
/* CRC and serial framing                                              */
/* ------------------------------------------------------------------ */

#define IRPROTO_SOF0 0xA5u
#define IRPROTO_SOF1 0x5Au
#define IRPROTO_OVERHEAD 7u      /* 5 header bytes + 2 CRC bytes */
#define IRPROTO_MAX_PAYLOAD 65535u

#define IRPROTO_ERR_SOF   (-10)
#define IRPROTO_ERR_LEN   (-11)
#define IRPROTO_ERR_CRC   (-12)
#define IRPROTO_ERR_SHORT (-13)
#define IRPROTO_ERR_CAP   (-14)

typedef struct {
    uint8_t type;
    uint16_t len;
    const uint8_t *payload;   /* points into the input buffer */
} irproto_frame_t;

/* CRC-16/CCITT-FALSE.  "123456789" -> 0x29B1. */
uint16_t ir_crc16_ccitt(const uint8_t *data, size_t n);

/* Packs one frame; returns IR_OK or a negative IRPROTO_ERR_*. */
int irproto_pack(uint8_t type, const uint8_t *payload, size_t len,
                 uint8_t *out, size_t cap, size_t *out_len);

/*
 * Parses one frame starting at buf[0].  On success *consumed is set to the
 * number of bytes the frame occupied.  Returns IR_OK or a negative error;
 * IRPROTO_ERR_SHORT means "not enough bytes yet", which is not corruption.
 */
int irproto_parse(const uint8_t *buf, size_t n, irproto_frame_t *out,
                  size_t *consumed);

#ifdef __cplusplus
}
#endif

#endif /* FRAME_H */
