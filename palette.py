"""Colors of the interface, the gel and the Excel export (spec §9).

Plain data without Qt, so that the Excel export and the command line can use
it too. Every widget takes its colors from here; nothing else hard-codes a color.
"""

LIGHT = dict(
    window='#F3F3F3', panel='#FFFFFF', line='#E0E0E0', text='#1F1F1F', text2='#616161',
    accent='#185FA5', on_accent='#FFFFFF', selected_bg='#E6F1FB', selected_text='#0C447C',
    robust_fg='#3B6D11', robust_bg='#EAF3DE', nominal_fg='#854F0B', nominal_bg='#FAEEDA',
    none_fg='#616161', none_bg='#EDEDED', error_fg='#A32D2D', error_bg='#FCEBEB',
    ref='#534AB7', alt='#993556', snp_bg='#FAC775', snp_text='#412402',
)
DARK = dict(
    window='#202020', panel='#2B2B2B', line='#3D3D3D', text='#F3F3F3', text2='#ABABAB',
    accent='#85B7EB', on_accent='#1F1F1F', selected_bg='#0C447C', selected_text='#E6F1FB',
    robust_fg='#C0DD97', robust_bg='#27500A', nominal_fg='#FAC775', nominal_bg='#633806',
    none_fg='#C8C8C8', none_bg='#3D3D3D', error_fg='#F7C1C1', error_bg='#791F1F',
    ref='#AFA9EC', alt='#ED93B1', snp_bg='#854F0B', snp_text='#FAEEDA',
)


def colors(dark: bool) -> dict:
    return dict(DARK if dark else LIGHT)


# The gel looks like a photo of a gel in both themes.
GEL = dict(background='#161616', well='#5B6670', band='#F2F2F2', marker='#8A8A8A', faint='#6A6A6A',
           label='#C8C8C8')

EXCEL_STATUS_FILLS = {'robust': 'EAF3DE', 'nominal': 'FAEEDA', 'none': 'EDEDED', 'error': 'FCEBEB'}
EXCEL_HEADER_FILL = 'EDEDED'
