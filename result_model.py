"""Russian titles of the core result columns and statuses."""
from i18n import N_

LABELS = {
    'assay_rank': N_('Ранг для SNP'), 'primer_pair_index': N_('Пара Primer3 №'),
    'primer_pair_penalty': N_('Штраф Primer3'), 'pairs_evaluated': N_('Проверено пар'),
    'genotype_quality': N_('Различимость генотипов'),
    'genotype_margin': N_('Оценка различимости (> 1)'),
    'worst_margin': N_('Оценка при ухудшении (> 1)'),
    'weakest_genotypes': N_('Труднее всего различить'),
    'bands_ref_ref': N_('Полосы REF/REF'), 'bands_ref_alt': N_('Полосы REF/ALT'),
    'bands_alt_alt': N_('Полосы ALT/ALT'),
    'worst_bands_ref_ref': N_('Полосы REF/REF при ухудшении'),
    'worst_bands_ref_alt': N_('Полосы REF/ALT при ухудшении'),
    'worst_bands_alt_alt': N_('Полосы ALT/ALT при ухудшении'),
    'visible_band_count': N_('Число полос в гетерозиготе'),
    'hidden_ref_count': N_('Невидимых фрагментов REF'),
    'hidden_alt_count': N_('Невидимых фрагментов ALT'),
    'search_mode': N_('Режим поиска'), 'scoring_model': N_('Модель оценки'),
    'variant': N_('Вариант'), 'mapped_id': N_('Контиг'), 'enzyme': N_('Фермент'),
    'site': N_('Сайт'), 'pattern': N_('Паттерн'), 'frags_ref': N_('Фрагменты REF, п.н.'),
    'frags_alt': N_('Фрагменты ALT, п.н.'), 'diag_delta_bp': N_('Различие, п.н.'),
    'amplicon_start': N_('Начало (1-based)'), 'amplicon_end': N_('Конец (1-based)'),
    'amplicon_len': N_('Длина, п.н.'), 'snp_offset_in_amplicon': 'SNP (0-based)',
    'primers_left': N_('Левый праймер 5′→3′'), 'primers_right': N_('Правый праймер 5′→3′'),
    'tm_left': N_('Tm лев., °C'), 'tm_right': N_('Tm прав., °C'), 'product_size': N_('Продукт, п.н.'),
    'primer3_size': N_('Продукт Primer3'), 'primer_left_start': N_('Начало L'),
    'primer_left_len': N_('Длина L'), 'primer_right_start': N_('Начало R'),
    'primer_right_len': N_('Длина R'), 'suppliers': N_('Поставщики'), 'status': N_('Статус'),
    'reason': N_('Пояснение'), 'delta_threshold_bp': N_('Порог Δ, п.н.'),
    'isoschizomers': N_('Изошизомеры'), 'primer_site': N_('Сайт в праймере'), 'repeat_pct': N_('Повторы, %'),
    'background_variants': N_('Известных вариантов в ампликоне'),
    'background_site_variants': N_('Варианты, меняющие разрезы/полосы'),
    'design': N_('Тип схемы'), 'primer_mismatches': N_('Несовпадения в праймере'),
    'line': N_('Строка'), 'chrom': N_('Хромосома'), 'pos': N_('Позиция'), 'ref': 'REF',
    'alt': 'ALT', 'details': N_('Подробности'), 'raw': N_('Исходная строка'),
}


STATUS_LABELS = {
    'no_discriminating_assay': 'Схема для трёх генотипов не найдена',
    'ok': 'Праймеры подобраны', 'window_only': 'Поиск по окну',
    'no_primers': 'Праймеры не найдены', 'primer3_error': 'Ошибка Primer3',
    'no_enzyme_found': 'Фермент не найден', 'analysis_error': 'Ошибка анализа',
}
