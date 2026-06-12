# PDL::Graphics::Cairo — API改善検討事項

## 問題1: subplots() の戻り値が行数で変わる (最優先)

### 現状
```perl
# rows=1 のとき: フラット戻り
my ($fig, $ax0, $ax1, ...) = subplots(1, N);

# rows>1 のとき: 行ごとの配列リファレンス
my ($fig, @rows) = subplots(M, N);  # M>1
$rows[$r][$c]    # Axesへのアクセス
```

### 問題点
- 同じ関数で戻り値の構造が変わるのは**驚き最小の原則**に違反
- `subplots(1,8)` が「横並び8列」なのか「縦積み8行」なのか直感に反する
- 作った本人でも間違える（EEGビューア開発中に2回以上ミス）

### 解決案A: 常にフラット戻り（推奨）
```perl
my ($fig, @ax) = subplots($nrows, $ncols, ...);
# インデックス: $ax[$r * $ncols + $c]
# ヘルパー: $fig->ax($r, $c) でアクセス
```

### 解決案B: 常に2D配列リファレンス
```perl
my ($fig, $axes) = subplots($nrows, $ncols, ...);
# $axes->[$r][$c]
```

### 解決案C: エイリアス `vstack(N)` / `hstack(N)` を追加
```perl
my ($fig, @ax) = vstack(8);  # 8行1列, フラット戻り
my ($fig, @ax) = hstack(8);  # 1行8列, フラット戻り
```

---

## 問題2: `add_subplot` が存在しない

matplotlib に慣れたユーザーが最初に試す `$fig->add_subplot($nrows, $ncols, $idx)` がない。
エラーメッセージも不親切。

### 解決案
```perl
# Figure.pm にエイリアス追加
sub add_subplot {
    my ($self, $nrows, $ncols, $idx) = @_;
    return $self->subplot($nrows, $ncols, $idx);
}
```

---

## 問題3: スライダーが上に振り切れる表示バグ

giza_server でスライダーのつまみ位置が上端に固定されて見える。
→ GSP_MSG_SLIDER / RESIZE の描画タイミング問題？ server13 以降の課題。

---

## 問題4: 複数ウィンドウがタブにまとまらない

`subplots(1, 8)` (1行8列) が8個の独立ウィンドウを開いた（最小化状態）。
→ `subplots` が内部で `subplot()` を8回呼び、各呼び出しで新規ウィンドウを
  開いてしまった可能性。giza_server のタブグループ化 (PID-based) の検証が必要。

