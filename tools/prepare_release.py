"""准备已确认的 v1.1.0 附件；只用 Python 标准库，不执行公开发布。"""

import argparse
import hashlib
import json
from pathlib import Path


# 故意绑定本次审阅后的版本和源码；后续发布仍须另行取得用户确认。
VERSION = "1.1.0"
SOURCE_SHA256 = "dbb3f527d39b7eab43c2e9c4d13e54adaddd3efc5a68e6103079b753ea68fe41"
ASSET_NAMES = ("focus_period.mpl", "快速使用.txt", "SHA256SUMS.txt")


def prepare(output: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    source = (root / "src/focus_period.mpl").read_bytes()
    if (root / "VERSION").read_text(encoding="utf-8").strip() != VERSION:
        raise ValueError("版本不属于本次已确认的 v1.1.0 发布")
    if hashlib.sha256(source).hexdigest() != SOURCE_SHA256:
        raise ValueError("源码与本次确认的发行字节不一致，停止打包")

    guide = (root / "docs/quick_start.txt").read_bytes()
    changelog = (root / "CHANGELOG.md").read_text(encoding="utf-8")
    details = changelog.split("## 1.1.0 — 2026-10-24\n\n", 1)[1]
    details = details.split('<a id="docs-20261020">', 1)[0].strip()
    quick_code = guide.decode("utf-8").split("\nrestart:\n", 1)[1]
    quick_code = "restart:\n" + quick_code.split("FocusPeriod:-Show(ans);", 1)[0]
    quick_code += "FocusPeriod:-Show(ans);\n"

    # 创建新目录，防止把已有发布材料静默覆盖。输出位置由调用方明确指定。
    output.mkdir(parents=True, exist_ok=False)
    (output / ASSET_NAMES[0]).write_bytes(source)
    (output / ASSET_NAMES[1]).write_bytes(guide)
    sums = "".join(
        f"{hashlib.sha256((output / name).read_bytes()).hexdigest()}  {name}\n"
        for name in ASSET_NAMES[:2]
    )
    (output / ASSET_NAMES[2]).write_bytes(sums.encode("utf-8"))
    notes = (
        "只运行计算程序时，下载附件 **focus_period.mpl** 即可；"
        "`快速使用.txt` 给出调用和查询方法，`SHA256SUMS.txt` 用于核对附件字节。"
        "需要全部示例、原生 `.mw` 工作表和说明时，下载 **Source code (zip)**。\n\n"
        "## 快速使用\n\n在 Maple 的可执行输入区运行，并修改 read 的实际路径：\n\n"
        f"```maple\n{quick_code}```\n\n"
        "**Windows 路径推荐使用 `/`；如用反斜杠，每个分隔符必须写两次。** "
        "先解决 `currentdir/read` 的第一个错误，再调用模块。"
        "红色的一维 Maple Input 可直接执行，无需转成二维公式；普通文字区不执行。"
        "替换旧文件后须 `restart` 并重新 `read`，确认版本为 `1.1.0`。\n\n"
        "最后两个数字分别为 g 与 p 的最高指数，可以独立设置。"
        "中间量可查看 `ans[\"G\"][2]`、`ans[\"H\"][0]`、"
        "`ans[\"u\"][3]`、`ans[\"g\"][3]`、`ans[\"p\"][2]`。\n\n"
        "[详细使用与变量表](https://github.com/Zerozero05/Maple_Focus_Period/blob/v1.1.0/README.md) · "
        "[高级用法](https://github.com/Zerozero05/Maple_Focus_Period/blob/v1.1.0/docs/advanced_usage.md) · "
        "[验证与限制](https://github.com/Zerozero05/Maple_Focus_Period/blob/v1.1.0/docs/validation.md)\n\n"
        f"## 本次更新、兼容性与验证\n\n{details}\n"
    )
    (output / "release-body.md").write_bytes(notes.encode("utf-8"))
    print(json.dumps({
        "version": VERSION,
        "assets": [{"name": name, "size": (output / name).stat().st_size,
                    "sha256": hashlib.sha256((output / name).read_bytes()).hexdigest()}
                   for name in ASSET_NAMES],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path, help="用于保存附件的新目录")
    prepare(parser.parse_args().output_dir)
