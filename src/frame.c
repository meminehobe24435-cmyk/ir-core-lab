/* frame.c -- BT.656 style timing and the serial frame protocol, see frame.h. */

#include "frame.h"

#include <string.h>

/* ------------------------------------------------------------------ */
/* Timing                                                              */
/* ------------------------------------------------------------------ */

void frame_preset_ir128(frame_timing_cfg_t *c)
{
    if (c == NULL) {
        return;
    }
    c->h_active = 128;
    c->h_blank  = 32;
    c->v_active = 128;
    c->v_blank  = 16;
    c->clocks_per_pixel = 1;
}

void frame_preset_bt656_625(frame_timing_cfg_t *c)
{
    if (c == NULL) {
        return;
    }
    c->h_active = 720;
    c->h_blank  = 144;
    c->v_active = 576;
    c->v_blank  = 49;
    c->clocks_per_pixel = 2;
}

long frame_line_clocks(const frame_timing_cfg_t *c)
{
    if (c == NULL || c->clocks_per_pixel <= 0) {
        return 0;
    }
    return (long)(c->h_active + c->h_blank) * (long)c->clocks_per_pixel;
}

long frame_total_clocks(const frame_timing_cfg_t *c)
{
    long lc = frame_line_clocks(c);
    if (c == NULL) {
        return 0;
    }
    return lc * (long)(c->v_active + c->v_blank);
}

void frame_timing_init(frame_timing_t *st, const frame_timing_cfg_t *c)
{
    if (st == NULL) {
        return;
    }
    memset(st, 0, sizeof(*st));
    st->cfg = c;
    if (c != NULL) {
        st->total_clocks = frame_total_clocks(c);
        st->line_clocks  = frame_line_clocks(c);
    }
}

int frame_timing_step(frame_timing_t *st, frame_signal_t *sig)
{
    long c, line_in_frame, clock_in_line;
    int line;

    if (st == NULL || sig == NULL || st->cfg == NULL) {
        return 0;
    }
    if (st->total_clocks <= 0 || st->line_clocks <= 0) {
        return 0;
    }

    c = st->frame_clocks_done;
    line_in_frame = c / st->line_clocks;
    clock_in_line = c % st->line_clocks;

    line = (int)line_in_frame;
    sig->clk = c;
    sig->y = line;
    sig->x = (int)(clock_in_line / (long)st->cfg->clocks_per_pixel);
    sig->de = (line < st->cfg->v_active && sig->x < st->cfg->h_active) ? 1 : 0;
    sig->active_pixel =
        (sig->de && (clock_in_line % (long)st->cfg->clocks_per_pixel == 0L)) ? 1 : 0;
    sig->hsync = (sig->x >= st->cfg->h_active) ? 1 : 0;
    sig->vsync = (line >= st->cfg->v_active) ? 1 : 0;

    if (sig->active_pixel) {
        st->active_pixels++;
    }
    if (clock_in_line == 0L) {
        st->lines_done++;
        if (line < st->cfg->v_active) {
            st->active_lines_done++;
        }
    }

    st->frame_clocks_done = c + 1;
    if (st->frame_clocks_done >= st->total_clocks) {
        st->last_active_pixels = st->active_pixels;
        st->last_lines         = st->lines_done;
        st->last_active_lines  = st->active_lines_done;
        st->frames_done++;
        st->frame_clocks_done  = 0L;
        st->active_pixels      = 0L;
        st->lines_done         = 0L;
        st->active_lines_done  = 0L;
        return 1;
    }
    return 0;
}

int frame_emit(const frame_timing_cfg_t *c, const uint8_t *img,
               int img_w, int img_h, uint8_t *stream, size_t cap,
               frame_report_t *rep)
{
    frame_timing_t st;
    frame_signal_t sig;

    if (c == NULL || c->h_active <= 0 || c->v_active <= 0 ||
        c->clocks_per_pixel <= 0) {
        return IR_ERR_PARAM;
    }
    if (rep != NULL) {
        memset(rep, 0, sizeof(*rep));
    }
    if (img != NULL && (img_w != c->h_active || img_h != c->v_active)) {
        /* Refuse rather than emit a raster that does not match the image. */
        return FRAME_ERR_DIM;
    }
    if (rep != NULL) {
        rep->dims_match = 1;
    }

    frame_timing_init(&st, c);
    if (stream != NULL && cap > 0u) {
        /* Blanking level first: every clock without data enable keeps 0x00 so
         * that the stream has no uninitialised holes. */
        memset(stream, 0, cap);
    }
    for (;;) {
        long slot = st.frame_clocks_done;
        int finished = frame_timing_step(&st, &sig);

        if (rep != NULL) {
            rep->clocks_emitted++;
            if (sig.de) {
                rep->active_pixel_clocks++;
            } else {
                rep->blanking_clocks++;
            }
            if (sig.hsync) {
                rep->hsync_asserted++;
            }
            if (sig.vsync) {
                rep->vsync_asserted++;
            }
        }
        if (sig.active_pixel && img != NULL) {
            uint8_t px = img[(size_t)sig.y * (size_t)img_w + (size_t)sig.x];
            if (stream != NULL && (size_t)slot < cap) {
                stream[slot] = px;
            }
            if (rep != NULL) {
                rep->pixels_emitted++;
            }
        }
        if (finished) {
            break;
        }
    }

    if (rep != NULL) {
        long active_clocks_expected =
            (long)c->h_active * (long)c->v_active * (long)c->clocks_per_pixel;
        int pixels_ok = (img == NULL) ? 1
                                      : (rep->pixels_emitted ==
                                         (long)c->h_active * (long)c->v_active);
        rep->total_clocks  = st.total_clocks;
        rep->lines         = st.last_lines;
        rep->active_lines  = st.last_active_lines;
        rep->counts_ok =
            (rep->clocks_emitted == st.total_clocks &&
             rep->active_pixel_clocks == active_clocks_expected &&
             rep->lines == (long)(c->v_active + c->v_blank) &&
             rep->active_lines == (long)c->v_active &&
             pixels_ok)
                ? 1
                : 0;
    }
    return IR_OK;
}

/* ------------------------------------------------------------------ */
/* CRC-16/CCITT-FALSE                                                  */
/* ------------------------------------------------------------------ */

uint16_t ir_crc16_ccitt(const uint8_t *data, size_t n)
{
    uint16_t crc = 0xFFFFu;
    size_t i;
    int b;

    if (data == NULL) {
        return crc;
    }
    for (i = 0u; i < n; ++i) {
        crc ^= (uint16_t)((uint16_t)data[i] << 8);
        for (b = 0; b < 8; ++b) {
            if (crc & 0x8000u) {
                crc = (uint16_t)((uint16_t)(crc << 1) ^ 0x1021u);
            } else {
                crc = (uint16_t)(crc << 1);
            }
        }
    }
    return crc;
}

/* ------------------------------------------------------------------ */
/* Serial framing                                                      */
/* ------------------------------------------------------------------ */

int irproto_pack(uint8_t type, const uint8_t *payload, size_t len,
                 uint8_t *out, size_t cap, size_t *out_len)
{
    size_t total;
    uint16_t crc;

    if (out == NULL || len > (size_t)IRPROTO_MAX_PAYLOAD) {
        return IR_ERR_PARAM;
    }
    if (payload == NULL && len > 0u) {
        return IR_ERR_PARAM;
    }
    total = len + IRPROTO_OVERHEAD;
    if (cap < total) {
        return IRPROTO_ERR_CAP;
    }

    out[0] = IRPROTO_SOF0;
    out[1] = IRPROTO_SOF1;
    out[2] = type;
    out[3] = (uint8_t)(len & 0xFFu);
    out[4] = (uint8_t)((len >> 8) & 0xFFu);
    if (len > 0u) {
        memcpy(&out[5], payload, len);
    }
    crc = ir_crc16_ccitt(&out[2], len + 3u);
    out[5 + len] = (uint8_t)(crc & 0xFFu);
    out[6 + len] = (uint8_t)((crc >> 8) & 0xFFu);

    if (out_len != NULL) {
        *out_len = total;
    }
    return IR_OK;
}

int irproto_parse(const uint8_t *buf, size_t n, irproto_frame_t *out,
                  size_t *consumed)
{
    size_t len, total;
    uint16_t crc_rx, crc_calc;

    if (buf == NULL || out == NULL) {
        return IR_ERR_PARAM;
    }
    if (n < 2u) {
        return IRPROTO_ERR_SHORT;
    }
    if (buf[0] != IRPROTO_SOF0 || buf[1] != IRPROTO_SOF1) {
        return IRPROTO_ERR_SOF;
    }
    if (n < 5u) {
        return IRPROTO_ERR_SHORT;
    }
    len = (size_t)buf[3] | ((size_t)buf[4] << 8);
    total = len + IRPROTO_OVERHEAD;
    if (n < total) {
        return IRPROTO_ERR_SHORT;
    }

    crc_rx = (uint16_t)buf[5 + len] | ((uint16_t)buf[6 + len] << 8);
    crc_calc = ir_crc16_ccitt(&buf[2], len + 3u);
    if (crc_rx != crc_calc) {
        return IRPROTO_ERR_CRC;
    }

    out->type = buf[2];
    out->len = (uint16_t)len;
    out->payload = &buf[5];
    if (consumed != NULL) {
        *consumed = total;
    }
    return IR_OK;
}
