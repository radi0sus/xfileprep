// sg-format.js
// One notation for space-group and point-group symbols everywhere on screen,
// as in International Tables: letters (lattice symbol, mirror and glide
// planes) in italics, digits upright, the screw component as a subscript
// (2₁, 4₂, 6₃) and the bar of a rotoinversion over the digit (1̄, 3̄, 4̄).
//
//   formatSG(hm)        HTML, for innerHTML     "P 1 21/c 1" → P 1 2₁/c 1
//   formatSGText(hm)    plain Unicode, for title="…", <option> and canvas
//   formatPG(symbol)    point groups / Laue classes ("2/m", "-3m", "222"):
//                       same style, but digits are never read as screw axes
//
// Only for display: the TITL line and all file contents keep the plain
// ASCII symbol with spaces (P 1 21/c 1).

const _SG_SUB = { '1': '₁', '2': '₂', '3': '₃', '4': '₄', '5': '₅' };
const _SG_SCREW = { '2': '1', '3': '12', '4': '123', '6': '12345' };

function _escSG(s) {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

// One space-separated token of an H-M symbol ("P", "21/a", "-4", "mm2" …).
// A two-digit run "nm" with m < n is a screw axis n_m (21, 31, 32, 41–43,
// 61–65); a point-group token such as "222" is formatted with screw = false.
function _sgToken(tok, html, screw) {
  let out = '';
  for (let i = 0; i < tok.length; i++) {
    const ch = tok[i];
    if (ch === '-' && /[1-6]/.test(tok[i + 1] || '')) {
      const d = tok[++i];
      out += html ? `<span class="sgBar">${d}</span>` : d + '̅';
    } else if (/[0-9]/.test(ch)) {
      const next = tok[i + 1] || '';
      if (screw && next && _SG_SCREW[ch] && _SG_SCREW[ch].includes(next) && !/[0-9]/.test(tok[i + 2] || '')) {
        out += html ? `${ch}<sub>${next}</sub>` : ch + _SG_SUB[next];
        i++;
      } else out += ch;
    } else if (/[A-Za-z]/.test(ch)) {
      out += html ? `<i>${ch}</i>` : ch;
    } else out += html ? _escSG(ch) : ch;
  }
  return out;
}

function _sgFormat(sym, html, screw) {
  if (sym === null || sym === undefined) return '';
  const s = String(sym).trim();
  // origin / axes suffix of table entries ("F d -3 m :2", "R -3 :H") stays upright
  const m = s.match(/^(.*?)\s*:\s*(\S+)$/);
  const body = m ? m[1] : s;
  const suffix = m ? `:${m[2]}` : '';
  // The symbol may come with or without spaces ("P 1 21/c 1", "P21/c");
  // unspaced symbols keep their compact form. A compact one ("P212121",
  // "P3121") is split with the Hall table, since digit runs are ambiguous
  // there (P312 has no screw axis, P3112 does).
  let tokens = body.split(/\s+/);
  const compact = screw && tokens.length === 1 && body.length > 2;
  if (compact && typeof HALL_SYMBOLS_TABLE !== 'undefined') {
    const hit = HALL_SYMBOLS_TABLE.find(e => e.hm.replace(/:.*$/, '').replace(/\s+/g, '') === body);
    if (hit) tokens = hit.hm.replace(/:.*$/, '').trim().split(/\s+/);
    else { const mm = body.match(/^([PABCFIRH])(.+)$/); if (mm) tokens = [mm[1], mm[2]]; }
  }
  const parts = tokens.map(t => _sgToken(t, html, screw));
  // thin spaces on screen (ITA prints the symbol without gaps), plain spaces in text
  const joined = parts.join(compact ? '' : (html ? '\u2009' : ' '));
  const text = joined + (suffix ? ` ${html ? _escSG(suffix) : suffix}` : '');
  return html ? `<span class="sgSym">${text}</span>` : text;
}

function formatSG(hm) { return _sgFormat(hm, true, true); }
function formatSGText(hm) { return _sgFormat(hm, false, true); }
function formatPG(symbol) { return _sgFormat(symbol, true, false); }
function formatPGText(symbol) { return _sgFormat(symbol, false, false); }

if (typeof module !== 'undefined') module.exports = { formatSG, formatSGText, formatPG, formatPGText };
