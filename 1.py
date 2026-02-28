import cairosvg
from PIL import Image
import io

# 读取 SVG 并转换为 PNG（带透明）
with open('BiomicsLOGO.svg', 'rb') as f:
    svg_data = f.read()

png_data = cairosvg.svg2png(bytestring=svg_data, output_width=500)

# 打开图像并添加白色背景
img = Image.open(io.BytesIO(png_data))

# 创建白色背景
background = Image.new('RGBA', img.size, (255, 255, 255, 255))

# 合并（如果原图有透明通道）
if img.mode == 'RGBA':
    background.paste(img, mask=img.split()[3])
else:
    background.paste(img)

# 保存为 RGB（去掉透明通道）
background.convert('RGB').save('BiomicsLOGO.png', 'PNG')