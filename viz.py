from PIL import Image,ImageDraw,ImageFont
import sys
import os 

def create_dir(directory):
    """"""
    path = os.getcwd().split('/')
    full_path = os.path.join('/'.join(path[:-2]), directory)
    os.mkdir(full_path)
    return full_path

# print(f'Directory created: {create_dir("hairpins_folder")}')

width_img = 800
height_img = 500
hairpin_string = "GGTGGACGAAGCCACCAFAFAASFASFASF"
stem_string = ["GGTGG",'CCACC']
loop_string = "ATGT"
# def hairpin_drawer(hairpin_string, stem_string, loop_string, width = 400, height = 300):
font = ImageFont.truetype("COURIER.ttf",30)
img  = Image.new(mode = "RGB", size = (width_img, height_img),color = 'white' )
draw = ImageDraw.Draw(img)

ascent, descent = font.getmetrics()
print(ascent,descent)

# Calculate the coordinates for centering the text
x = (img.width) // 2 - ascent
y = (img.height) // 2 - descent

draw.text((x/2,y), hairpin_string, fill='Blue',font = font,antialias=True, align='centre')

img.show()
# img.save(f"hairpins_folder/image_text.jpg")


